import os
from typing import Dict, List
import urllib.request
from itertools import cycle

from mmcif.io.IoAdapterCore import IoAdapterCore  # type: ignore
from mmcif.api.PdbxContainers import DataContainer  # type: ignore

import molviewspec as mvs

mmcif_url = lambda x: f"https://files.rcsb.org/download/{x}.cif"


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
    The default view of a structure from the pdb.

    Each polymer entity is represented as cartoon.
    Each branched entity is represented as a carbohydrate.
    All other entities are represented by  ball_and_stick.
    Colors are different for each entity.
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
