import os
from typing import Dict, List, Literal
import urllib.request
from itertools import cycle

from mmcif.io.IoAdapterCore import IoAdapterCore  # type: ignore
from mmcif.api.PdbxContainers import DataContainer  # type: ignore

import molviewspec as mvs
from molviewspec import Snapshot

mmcif_url = lambda x: f"https://files.rcsb.org/download/{x}.cif"


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


# TODO: use streaming IoAdapterPy and cache
def get_data_container(pdb_id: str) -> DataContainer:
    """Download and parse an mmCIF file for a given PDB ID.

    Args:
        pdb_id: PDB identifier (e.g., '1abc')

    Returns:
        Parsed mmCIF data container
    """

    url = mmcif_url(pdb_id)
    filepath_local = os.path.join("./files", f"{pdb_id}.cif")
    urllib.request.urlretrieve(url, filepath_local)
    io = IoAdapterCore()
    data_container = io.readFile(filepath_local)[0]
    return data_container


def get_entity_types_by_id(pdb_id: str) -> Dict[str, str]:
    """Get entity types mapped by entity ID for a PDB structure.

    Args:
        pdb_id: PDB identifier

    Returns:
        Dictionary mapping entity ID to entity type (polymer, water, etc.)
    """
    data_container = get_data_container(pdb_id)
    entities = data_container.getObj("entity")
    entities.getAttributeValueList("type")
    type_by_id = dict(
        zip(
            entities.getAttributeValueList("id"), entities.getAttributeValueList("type")
        )
    )
    return type_by_id


def get_sequence_metadata(pdb_id: str):
    """Extract combined polymer and non-polymer chain metadata including sequence numbering.

    Args:
        pdb_id: PDB identifier

    Returns:
        Dictionary with chain IDs as keys and metadata (auth_seq_num, pdb_seq_num, pdb_mon_id) as values.
        Combines data from both polymer and non-polymer chains.
    """
    dc = get_data_container(pdb_id)
    combined_data = {}

    # Get polymer data (uses pdb_strand_id)
    try:
        poly_scheme = dc.getObj("pdbx_poly_seq_scheme")
        poly_chain_ids = poly_scheme.getAttributeUniqueValueList("pdb_strand_id")
        for chain_id in poly_chain_ids:
            pdb_seq_nums = poly_scheme.selectValuesWhereConditions(
                attributeName="pdb_seq_num", conditionsD={"pdb_strand_id": chain_id}
            )
            auth_seq_nums = poly_scheme.selectValuesWhereConditions(
                attributeName="auth_seq_num", conditionsD={"pdb_strand_id": chain_id}
            )
            pdb_mon_ids = poly_scheme.selectValuesWhereConditions(
                attributeName="pdb_mon_id", conditionsD={"pdb_strand_id": chain_id}
            )
            combined_data[chain_id] = {
                "auth_seq_num": [int(i) for i in auth_seq_nums],
                "pdb_seq_num": [int(i) for i in pdb_seq_nums],
                "pdb_mon_id": pdb_mon_ids,
            }
    except:
        # No polymer data available
        pass

    # Get non-polymer data (uses asym_id)
    try:
        nonpoly_scheme = dc.getObj("pdbx_nonpoly_scheme")
        nonpoly_chain_ids = nonpoly_scheme.getAttributeUniqueValueList("asym_id")
        for chain_id in nonpoly_chain_ids:
            pdb_seq_nums = nonpoly_scheme.selectValuesWhereConditions(
                attributeName="pdb_seq_num", conditionsD={"asym_id": chain_id}
            )
            auth_seq_nums = nonpoly_scheme.selectValuesWhereConditions(
                attributeName="auth_seq_num", conditionsD={"asym_id": chain_id}
            )
            pdb_mon_ids = nonpoly_scheme.selectValuesWhereConditions(
                attributeName="pdb_mon_id", conditionsD={"asym_id": chain_id}
            )
            combined_data[chain_id] = {
                "auth_seq_num": [int(i) for i in auth_seq_nums],
                "pdb_seq_num": [int(i) for i in pdb_seq_nums],
                "pdb_mon_id": pdb_mon_ids,
            }
    except:
        # No non-polymer data available
        pass

    return combined_data


def get_residue_name(residue_selection: Dict[str, str | int], pdb_id: str) -> str:
    """Get the residue name (pdb_mon_id) for a given residue selection.

    Args:
        residue_selection: Dictionary with label_asym_id and label_seq_id
        pdb_id: PDB identifier

    Returns:
        Residue name (pdb_mon_id) as string

    Raises:
        ValueError: If residue selection is invalid or residue not found
    """
    metadata = get_sequence_metadata(pdb_id)

    chain_id = residue_selection["label_asym_id"]
    seq_id = residue_selection["label_seq_id"]

    if chain_id not in metadata:
        raise ValueError(f"Chain {chain_id} not found in structure")

    chain_metadata = metadata[chain_id]

    # Find the index of the sequence ID
    try:
        seq_index = chain_metadata["pdb_seq_num"].index(seq_id)
        return chain_metadata["pdb_mon_id"][seq_index]
    except ValueError:
        raise ValueError(f"Residue {seq_id} not found in chain {chain_id}")


def validate_single_residue_selection(selection: Dict[str, str | int], pdb_id: str):
    """Validate that a single residue selection is valid for the given structure.

    Args:
        selection: Dictionary with label_asym_id and label_seq_id
        pdb_id: PDB identifier

    Returns:
        True if selection is valid

    Raises:
        ValueError: If selection is invalid
    """
    metadata = get_sequence_metadata(pdb_id=pdb_id)

    chain_id = selection["label_asym_id"]
    seq_id = selection["label_seq_id"]
    selection_txt = f"{chain_id}:{seq_id}"

    # Check if chain exists
    if chain_id not in metadata:
        raise ValueError(
            f'Selection {selection_txt} not valid. Chain {chain_id} not available. Chains are: {",".join(metadata.keys())}'
        )

    chain_metadata = metadata[chain_id]

    # Check if sequence ID exists
    pdb_seq_id_ok = seq_id in chain_metadata["pdb_seq_num"]
    auth_seq_id_ok = seq_id in chain_metadata["auth_seq_num"]

    if (not pdb_seq_id_ok) and auth_seq_id_ok:
        raise ValueError(
            f"Residue selection {selection_txt} not valid with PDB numbering but valid with author numbering. Residue selections should use PDB numbering (zero indexed). Chain {chain_id} has {len(chain_metadata['pdb_seq_num'])} residues."
        )

    if not pdb_seq_id_ok:
        raise ValueError(
            f"Residue selection {selection_txt} not valid with PDB numbering or author numbering. Residue selections should use PDB numbering (zero indexed). Chain {chain_id} has {len(chain_metadata['pdb_seq_num'])} residues."
        )

    return True


def validate_polymer_range_selection(
    selection: Dict[str, str | int | None],
    metadata: Dict[str, Dict[str, List[str | int]]],
):
    """Validate that a residue range selection is valid for the given structure.

    Args:
        selection: Dictionary with label_asym_id, beg_label_seq_id, end_label_seq_id
        metadata: Structure metadata from get_sequence_metadata

    Returns:
        True if selection is valid

    Raises:
        ValueError: If selection is invalid
    """
    chain_id = selection["label_asym_id"]
    beg_id = selection["beg_label_seq_id"]
    end_id = selection["end_label_seq_id"]
    selection_txt = f"{chain_id}:{beg_id}-{end_id}"
    if not chain_id in metadata.keys():
        raise ValueError(
            f'Selection {selection_txt} not valid. Chain {chain_id} not available. Chains are: {",".join(metadata.keys())}'
        )

    pdb_seq_id_ok = (beg_id in metadata[chain_id]["pdb_seq_num"]) and (
        end_id in metadata[chain_id]["pdb_seq_num"]
    )

    auth_seq_id_ok = (beg_id in metadata[chain_id]["auth_seq_num"]) and (
        end_id in metadata[chain_id]["auth_seq_num"]
    )

    if (not pdb_seq_id_ok) and auth_seq_id_ok:
        raise ValueError(
            f"Residue selection {selection_txt} not valid with PDB numbering but valid with author numbering. Residue selections should use PDB numbering (zero indexed). Chain {chain_id} has {len(metadata[chain_id]['pdb_seq_num'])} residues."
        )

    if not pdb_seq_id_ok:
        raise ValueError(
            f"Residue selection {selection_txt} not valid with PDB numbering or author numbering. Residue selections should use PDB numbering (zero indexed). Chain {chain_id} has {len(metadata[chain_id]['pdb_seq_num'])} residues."
        )

    return True


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


def story_wrapper(snapshots: List[Snapshot], description: str | None = None):
    """Wrap snapshots into a story (States object) with metadata.

    Args:
        snapshots: List of snapshots to include in the story
        description: Optional description for the story

    Returns:
        JSON string of the States object
    """
    if not description:
        description = f"{len(snapshots)} story"

    states_json = mvs.States(
        snapshots=snapshots, metadata=mvs.GlobalMetadata(description=description)
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


def highlight_residues(
    pdb_id: str,
    selections_to_highlight: List[Dict[str, str | int | None]],
    title: str | None = None,
    description: str | None = None,
    highlight_representation: Literal["ball_and_stick", "cartoon"] = "cartoon",
    baseline_opacity: float = 0.9,
):
    """
    Creates a snapshot which highlights different sections of the polymer chains on a given pdb. The whole polymer is shown
    translucent, while the highlights are fully opaque. Only the asymmetric unit is shown due to an issue with label positioning.

    Args:
        pdb_id: PDB identifier
        selections_to_highlight: List of dictionaries defining residue selections.
            Each dictionary should have the structure:
            {
                "label_asym_id": str,  # Chain identifier (e.g., "A", "B")
                "beg_label_seq_id": int,  # Starting residue number (PDB numbering)
                "end_label_seq_id": int,  # Ending residue number (PDB numbering)
                "text_label": str | None,  # Optional text label for the selection
                "color": str | None,  # Optional color (e.g., "yellow", "red")
            }
        title: Optional title for the snapshot
        description: Optional description for the snapshot
        highlight_representation: Representation type for highlights
        baseline_opacity: Opacity for the remaining structure.
    """
    # Sort by selection length (smallest first)
    selections_sorted = sorted(
        selections_to_highlight,
        key=lambda x: x["end_label_seq_id"] - x["beg_label_seq_id"],
    )
    builder = mvs.create_builder()

    structure = (
        builder.download(url=mmcif_url(pdb_id)).parse(format="mmcif").model_structure()
    )
    type_by_id = get_entity_types_by_id(pdb_id)
    poly_by_id = {k: v for k, v in type_by_id.items() if v == "polymer"}

    poly_meta = get_sequence_metadata(pdb_id)

    # Add smallest selections first - these will be rendered first.
    for selection in selections_sorted:
        color = selection.get("color") or next(HIGHLIGHT_COLORS)
        if validate_polymer_range_selection(selection, poly_meta):
            (
                structure.component(selector=selection)
                .representation(type=highlight_representation)
                .color(color=color)
                # .opacity(opacity=highlight_opacity)
            )
        label_text = (
            selection.get("text_label")
            or f"{selection['beg_label_seq_id']}-{selection['end_label_seq_id']}: {selection['label_asym_id']}"
        )
        label_seq_id = int(
            (selection["beg_label_seq_id"] + selection["end_label_seq_id"]) / 2
        )
        label_selector = {
            "label_seq_id": label_seq_id,
            "label_asym_id": selection["label_asym_id"],
        }
        structure.component(selector=label_selector).label(text=label_text)

    # Add representation for whole protein
    for entity_id, entity_type in poly_by_id.items():
        color = next(COLORS)
        (
            structure.component(selector={"label_entity_id": entity_id})
            .representation(type=REP_BY_ENTITY_TYPE.get(entity_type, DEFAULT_REP))
            .color(color=color)
            .opacity(opacity=baseline_opacity)
        )

    snapshot = builder.get_snapshot(
        title=title or f"{pdb_id.upper()}",
        description=description
        or f"Highlight of {len(selections_to_highlight)} sections",
        description_format="markdown",
        linger_duration_ms=DEFAULT_LINGER_MS,
        transition_duration_ms=DEFAULT_TRANSITION_MS,
    )
    return snapshot


def residue_with_interactions(
    pdb_id: str,
    residue_selection: Dict[str, int | str],
    baseline_opacity: float = 0.5,
    title: str | None = None,
    description: str | None = None,
):
    """Create snapshot focusing on a residue with non-covalent interactions.

    Args:
        pdb_id: PDB identifier
        residue_selection: Dictionary with label_asym_id, label_seq_id, optional text_label and color
        baseline_opacity: Opacity for the background structure
        title: Optional snapshot title
        description: Optional snapshot description

    Returns:
        Snapshot object with residue interactions highlighted
    """
    type_by_id = get_entity_types_by_id(pdb_id)

    builder = mvs.create_builder()

    structure = (
        builder.download(url=mmcif_url(pdb_id))
        .parse(format="mmcif")
        .model_structure()
    )

    for entity_id, entity_type in type_by_id.items():
        if entity_type == "water":
            continue
        (
            structure.component(selector={"label_entity_id": entity_id})
            .representation(type=REP_BY_ENTITY_TYPE.get(entity_type, DEFAULT_REP))
            .color(color=next(COLORS))
            .opacity(opacity=baseline_opacity)
        )

    residue_selector = {
        "label_asym_id": residue_selection["label_asym_id"],
        "label_seq_id": residue_selection["label_seq_id"],
    }
    if validate_single_residue_selection(residue_selector, pdb_id):
        res_name = get_residue_name(residue_selector, pdb_id)

        color = residue_selection.get("color") or next(HIGHLIGHT_COLORS)
        selection_label = f"{res_name} {residue_selection['label_seq_id']} {residue_selection['label_asym_id']}"
        label_text = residue_selection.get("text_label") or selection_label
 
        (
            structure.component(
                selector=residue_selector,
                custom={
                    "molstar_show_non_covalent_interactions": True,
                    "molstar_non_covalent_interactions_radius_ang": 5.0,
                },
            )
            .representation(type="ball_and_stick")
            .color(color=color)
        )
        structure.component(selector=residue_selector).label(text=label_text).focus()

        snapshot = builder.get_snapshot(
            title=title or f"{pdb_id.upper()}",
            description=description or f"Residue interactions for {selection_label}",
            description_format="markdown",
            linger_duration_ms=DEFAULT_LINGER_MS,
            transition_duration_ms=DEFAULT_TRANSITION_MS,
        )
    return snapshot
