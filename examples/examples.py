import os
from typing import Dict, List
import urllib.request
from itertools import cycle

from mmcif.io.IoAdapterCore import IoAdapterCore  # type: ignore
from mmcif.api.PdbxContainers import DataContainer  # type: ignore

import molviewspec as mvs

mmcif_url = lambda x: f"https://files.rcsb.org/download/{x}.cif"


# TODO: use streaming IoAdapterPy and cache
def get_data_container(pdb_id: str) -> DataContainer:

    url = mmcif_url(pdb_id)
    filepath_local = os.path.join("./files", f"{pdb_id}.cif")
    urllib.request.urlretrieve(url, filepath_local)
    io = IoAdapterCore()
    data_container = io.readFile(filepath_local)[0]
    return data_container


def get_entity_types_by_id(pdb_id: str) -> Dict[str, str]:
    data_container = get_data_container(pdb_id)
    entities = data_container.getObj("entity")
    entities.getAttributeValueList("type")
    type_by_id = dict(
        zip(
            entities.getAttributeValueList("id"), entities.getAttributeValueList("type")
        )
    )
    return type_by_id


def get_polymer_metadata(pdb_id: str):
    dc = get_data_container(pdb_id)
    scheme = dc.getObj("pdbx_poly_seq_scheme")
    chain_ids = scheme.getAttributeUniqueValueList("pdb_strand_id")
    pdb_seq_nums = [
        scheme.selectValuesWhereConditions(
            attributeName="pdb_seq_num", conditionsD={"pdb_strand_id": chain_id}
        )
        for chain_id in chain_ids
    ]
    pdb_seq_nums = [[int(i) for i in seq] for seq in pdb_seq_nums]
    auth_seq_nums = [
        scheme.selectValuesWhereConditions(
            attributeName="auth_seq_num", conditionsD={"pdb_strand_id": chain_id}
        )
        for chain_id in chain_ids
    ]
    auth_seq_nums = [[int(i) for i in seq] for seq in auth_seq_nums]
    pdb_mon_id = [
        scheme.selectValuesWhereConditions(
            attributeName="pdb_mon_id", conditionsD={"pdb_strand_id": chain_id}
        )
        for chain_id in chain_ids
    ]
    data = {
        chain_id: {
            "auth_seq_num": auth_seq_nums[i],
            "pdb_seq_nums": pdb_seq_nums[i],
            "pdb_mon_id": pdb_mon_id[i],
        }
        for i, chain_id in enumerate(chain_ids)
    }

    return data


COLORS = cycle(
    [
        "blue",
        "darkgreen",
        "lightblue",
        "lightgreen",
        "palegoldenrod",
        "lightcoral",
        "purple",
    ]
)

HIGHLIGHT_COLORS = cycle(["yellow", "whitesmoke", "black"])

WATER_COLOR = "aqua"

REP_BY_ENTITY_TYPE = {
    "polymer": "cartoon",
    "non-polymer": "ball_and_stick",
    "branched": "carbohydrate",
    "water": "ball_and_stick",
    "macrolide": "ball_and_stick",
}

DEFAULT_REP = "ball_and_stick"
DEFAULT_LINGER_MS = 4000
DEFAULT_TRANSITION_MS = 2000


def default_story(
    pdb_ids: List[str],
    titles: List[str] | None = None,
    descriptions_md: List[str] | None = None,
    global_description: str | None = None,
) -> str:
    """
    Creates a default story story from a list of pdbs.
    """
    n_pdbs = len(pdb_ids)
    titles = titles or [None] * n_pdbs
    descriptions_md = descriptions_md or [None] * n_pdbs
    global_description = global_description or f'Snapshots of {", ".join(pdb_ids)}'
    snapshots = [
        default_snapshot(pdb_id, title, description)
        for pdb_id, title, description in zip(pdb_ids, titles, descriptions_md)
    ]
    states_json = mvs.States(
        snapshots=snapshots, metadata=mvs.GlobalMetadata(description=global_description)
    ).dumps(indent=2)
    return states_json


def default_snapshot(
    pdb_id: str, title: str | None = None, description: str | None = None
):
    """
    Creates a deafult snapshot for including in an animation (series of states)

    Each polymer entity is represented as cartoon.
    Each branched entity is represented as a carbohydrate.
    All other entities are represented by  ball_and_stick.
    Colors are different for each entity except water which is always the same color.
    """
    type_by_id = get_entity_types_by_id(pdb_id)

    builder = mvs.create_builder()

    structure = (
        builder.download(url=mmcif_url(pdb_id))
        .parse(format="mmcif")
        .assembly_structure()
    )

    for entity_id, entity_type in type_by_id.items():
        if entity_type == "water":
            color = WATER_COLOR
        else:
            color = next(COLORS)
        (
            structure.component(selector={"label_entity_id": entity_id})
            .representation(type=REP_BY_ENTITY_TYPE.get(entity_type, DEFAULT_REP))
            .color(color=color)
        )
    snapshot = builder.get_snapshot(
        title=title or f"{pdb_id.upper()}",
        description=description or "Default view",
        description_format="markdown",
        linger_duration_ms=DEFAULT_LINGER_MS,
        transition_duration_ms=DEFAULT_TRANSITION_MS,
    )
    return snapshot


class ResiduesSelection:
    text_label: str
    label_asym_id: str
    beg_label_seq_id: int
    end_label_seq_id: int

    def __lt__(self, obj):
        len_this = self.end_label_seq_id - self.beg_label_seq_id
        len_that = obj.end_label_seq_id - obj.beg_label_seq_id
        return len_this < len_that

def validate_selection(selection: ResiduesSelection, metadata: Dict[str, Dict[str, List[str | int]]]): 
    chain_id = selection.label_asym_id
    if not chain_id in metadata.keys(): 
        raise ValueError(f'Requested chain ID (label_asym_id): {selection.label_asym_id} not available. Available chains are: {",".join(metadata.keys())}')
    
    pdb_seq_id_ok = (selection.beg_label_seq_id in metadata[chain_id]['pdb_seq_num']) and (selection.end_label_seq_id in metadata[chain_id]['pdb_seq_num'])

    auth_seq_id_ok = (selection.beg_label_seq_id in metadata[chain_id]['auth_seq_num']) and (selection.end_label_seq_id in metadata[chain_id]['auth_seq_num'])

    if (not pdb_seq_id_ok) and (not auth_seq_id_ok): 
        raise ValueError('Requestion')
     


def highlight_residues(
    pdb_id: str,
    selections_to_hightlight: List[ResiduesSelection],
    title: str | None = None,
    description: str | None = None,
):
    """
    Creates a snapshot which highlights different sections of the polymer chains
    """
    selections_sorted = sorted(selections_to_hightlight)
    builder = mvs.create_builder()

    structure = (
        builder.download(url=mmcif_url(pdb_id))
        .parse(format="mmcif")
        .assembly_structure()
    )
    type_by_id = get_entity_types_by_id(pdb_id)
    poly_by_id = {k: v for k, v in type_by_id.items() if v == "polymer"}

    # Add smallest selections first
    for selections in selections_sorted:
        
    # Add representation for whole protein
    for entity_id, entity_type in poly_by_id.items():
        color = next(COLORS)
        (
            structure.component(selector={"label_entity_id": entity_id})
            .representation(type=REP_BY_ENTITY_TYPE.get(entity_type, DEFAULT_REP))
            .color(color=color)
        )

