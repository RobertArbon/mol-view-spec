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
    global_description: str,
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
            default_snapshot(
                pdb_id,
                pdb_config["type_by_id"],
                pdb_config["title"],
                pdb_config["description"],
            )
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
    description = f"# {title}\n{description}"
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


if __name__ == "__main__":
    # Example usage of default_story
    config = {
        "1cbs": {
            "title": "Crambin Structure",
            "description": "Small hydrophobic protein from *Crambe abyssinica*",
            "type_by_id": {"1": "polymer", "2": "water"},
        },
        "1ubq": {
            "title": "Ubiquitin Structure",
            "description": "Small regulatory protein found in eukaryotic cells",
            "type_by_id": {"1": "polymer", "2": "water"},
        },
    }

    # Create story using default_story
    story_json = default_story(config, "Comparison of Crambin and Ubiquitin structures")
