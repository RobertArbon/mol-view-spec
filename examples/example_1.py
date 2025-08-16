from typing import Dict, List, Literal
from itertools import cycle

import molviewspec as mvs
from molviewspec import Snapshot



def mmcif_url(pdb_id):
    return f"https://files.rcsb.org/download/{pdb_id}.cif"


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

PdbType = str
EntityIdType = str
EntityType = Literal["polymer", "non-polymer", "water", "branched"]

def default_story(
    config: Dict[PdbType, Dict[str, str | Dict[EntityIdType, EntityType]]], 
    global_description: str 
) -> str:
    """
    Creates a default story story from a list of pdbs.
    config = {
        pdb_id: str # pdb id
        title: title of the pdb
        description: markdown formatted description or other text accompanying this pdb.
        type_by_id: the entity type by the label_entity_id
    }
    """
    snapshots = []
    for pdb_id, pdb_config in config.items():  
        snapshots.append(
            default_snapshot(pdb_id, pdb_config['type_by_id'], pdb_config['title'], pdb_config['description'])
            )
    
    states_json = mvs.States(
        snapshots=snapshots, metadata=mvs.GlobalMetadata(description=global_description)
    ).dumps(indent=2)
    return states_json


def default_snapshot(
    pdb_id: str,
    type_by_id: Dict[str, EntityType],
    title: str,
    description: str,
):
    """
    Creates a deafult snapshot for including in an animation (series of states)

    Each polymer entity is represented as cartoon.
    Each branched entity is represented as a carbohydrate.
    All other entities are represented by  ball_and_stick.
    Colors are different for each entity except water which is always the same color.
    """
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
        title=title,  
        description=description, 
        description_format="markdown",
        linger_duration_ms=DEFAULT_LINGER_MS,
        transition_duration_ms=DEFAULT_TRANSITION_MS,
    )
    return snapshot


def story_wrapper(snapshots: List[Snapshot], description: str):
    """Wrap snapshots into a story (States object) with metadata.

    Args:
        snapshots: List of snapshots to include in the story
        description: Optional description for the story

    Returns:
        JSON string of the States object
    """

    states_json = mvs.States(
        snapshots=snapshots, metadata=mvs.GlobalMetadata(description=description)
    ).dumps(indent=2)
    return states_json


# def highlight_residues(
#     pdb_id: str,
#     selections_to_highlight: List[Dict[str, str | int | None]],
#     title: str | None = None,
#     description: str | None = None,
#     highlight_representation: Literal["ball_and_stick", "cartoon"] = "cartoon",
#     baseline_opacity: float = 0.9,
# ):
#     """
#     Creates a snapshot which highlights different sections of the polymer chains on a given pdb. The whole polymer is shown
#     translucent, while the highlights are fully opaque. Only the asymmetric unit is shown due to an issue with label positioning.

#     Args:
#         pdb_id: PDB identifier
#         selections_to_highlight: List of dictionaries defining residue selections.
#             Each dictionary should have the structure:
#             {
#                 "label_asym_id": str,  # Chain identifier (e.g., "A", "B")
#                 "beg_label_seq_id": int,  # Starting residue number (PDB numbering)
#                 "end_label_seq_id": int,  # Ending residue number (PDB numbering). Can be the same as beg_label_seq_id
#                 "text_label": str | None,  # Optional text label for the selection
#                 "color": str | None,  # Optional color (e.g., "yellow", "red")
#             }
#         title: Optional title for the snapshot
#         description: Optional description for the snapshot
#         highlight_representation: Representation type for highlights
#         baseline_opacity: Opacity for the remaining structure.
#     """
#     # Sort by selection length (smallest first)
#     selections_sorted = sorted(
#         selections_to_highlight,
#         key=lambda x: x["end_label_seq_id"] - x["beg_label_seq_id"],
#     )
#     builder = mvs.create_builder()

#     structure = (
#         builder.download(url=mmcif_url(pdb_id)).parse(format="mmcif").model_structure()
#     )
#     # Add smallest selections first - these will be rendered first.
#     for selection in selections_sorted:
#         color = selection.get("color") or next(HIGHLIGHT_COLORS)
#         (
#             structure.component(selector=selection)
#             .representation(type=highlight_representation)
#             .color(color=color)
#         )
#         label_text = (
#             selection.get("text_label")
#             or f"{selection['beg_label_seq_id']}-{selection['end_label_seq_id']}: {selection['label_asym_id']}"
#         )
#         label_seq_id = int(
#             (selection["beg_label_seq_id"] + selection["end_label_seq_id"]) / 2
#         )
#         label_selector = {
#             "label_seq_id": label_seq_id,
#             "label_asym_id": selection["label_asym_id"],
#         }
#         structure.component(selector=label_selector).label(text=label_text)

#     # Add representation for whole protein - showing all polymers
#     # Note: Without validation, we assume all entities are polymers
#     structure_component = structure.component(selector={})
#     for _ in range(1):  # Single iteration for basic representation
#         entity_type = "polymer"
#         color = next(COLORS)
#         (
#             structure_component
#             .representation(type=REP_BY_ENTITY_TYPE.get(entity_type, DEFAULT_REP))
#             .color(color=color)
#             .opacity(opacity=baseline_opacity)
#         )

#     snapshot = builder.get_snapshot(
#         title=title or f"{pdb_id.upper()}",
#         description=description
#         or f"Highlight of {len(selections_to_highlight)} sections",
#         description_format="markdown",
#         linger_duration_ms=DEFAULT_LINGER_MS,
#         transition_duration_ms=DEFAULT_TRANSITION_MS,
#     )
#     return snapshot


# def residue_with_interactions(
#     pdb_id: str,
#     residue_selection: Dict[str, int | str],
#     baseline_opacity: float = 0.5,
#     title: str | None = None,
#     description: str | None = None,
# ):
#     """Create snapshot focusing on a residue with non-covalent interactions.

#     Args:
#         pdb_id: PDB identifier
#         residue_selection: dictionary representing the selection criteria for a single residue:
#             {
#                 "label_asym_id": str,  # Chain identifier (e.g., "A", "B")
#                 "label_seq_id": int, # Residue number (pdb numbering)
#                 "text_label": str | None,  # Optional text label for the selection
#                 "color": str | None,  # Optional color (e.g., "yellow", "red")
#             }
#         baseline_opacity: Opacity for the background structure
#         title: Optional snapshot title
#         description: Optional snapshot description

#     Returns:
#         Snapshot object with residue interactions highlighted
#     """
#     builder = mvs.create_builder()

#     structure = (
#         builder.download(url=mmcif_url(pdb_id)).parse(format="mmcif").model_structure()
#     )

#     # Add basic structure representation without validation
#     (
#         structure.component(selector={})
#         .representation(type="cartoon")
#         .color(color=next(COLORS))
#         .opacity(opacity=baseline_opacity)
#     )

#     residue_selector = {
#         "label_asym_id": residue_selection["label_asym_id"],
#         "label_seq_id": residue_selection["label_seq_id"],
#     }
    
#     color = residue_selection.get("color") or next(HIGHLIGHT_COLORS)
#     selection_label = f"Residue {residue_selection['label_seq_id']} {residue_selection['label_asym_id']}"
#     label_text = residue_selection.get("text_label") or selection_label

#     (
#         structure.component(
#             selector=residue_selector,
#             custom={
#                 "molstar_show_non_covalent_interactions": True,
#                 "molstar_non_covalent_interactions_radius_ang": 5.0,
#             },
#         )
#         .representation(type="ball_and_stick")
#         .color(color=color)
#     )
#     structure.component(selector=residue_selector).label(text=label_text).focus()

#     snapshot = builder.get_snapshot(
#         title=title or f"{pdb_id.upper()}",
#         description=description or f"Residue interactions for {selection_label}",
#         description_format="markdown",
#         linger_duration_ms=DEFAULT_LINGER_MS,
#         transition_duration_ms=DEFAULT_TRANSITION_MS,
#     )
#     return snapshot

if __name__=='__main__': 
