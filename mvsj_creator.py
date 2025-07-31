import molviewspec as mvs


def create_animation():
    """
    Create a 2 state animation using the pdb 5vpy
    State 1:
    - title: 'State 1'
    - description: 'The protein and ligand in canonical representation'
    - show the protein in cartoon representation.
    - show the ligand in ball and stick representation
    State 2:
    - title: 'State 2'
    - description: 'The switch II region'
    - highlight the residues 22 - 44 in a different colour.
    """
    # State 1
    builder1 = mvs.create_builder()
    structure1 = builder1.download(url="https://files.rcsb.org/download/5vpy.cif").parse(format="mmcif").model_structure()
    structure1.component(selector="protein").representation(type="cartoon")
    structure1.component(selector="ligand").representation(type="ball_and_stick")
    state1 = builder1.get_state(title='State 1', description='The protein and ligand in canonical representation')

    # State 2
    builder2 = mvs.create_builder()
    structure2 = builder2.download(url="https://files.rcsb.org/download/5vpy.cif").parse(format="mmcif").model_structure()
    structure2.component(selector="protein").representation(type="cartoon")
    structure2.component(selector="ligand").representation(type="ball_and_stick")

    # In PDB entry 5vpy, the protein chain is A. We select residues 22-44 for highlighting.
    switch_ii_selector = {"label_asym_id": "A", "beg_label_seq_id": 22, "end_label_seq_id": 44}
    structure2.component(selector=switch_ii_selector).representation(type="cartoon", color="#FF0000")  # red highlight

    state2 = builder2.get_state(title='State 2', description='The switch II region')

    # Combine states into a multi-state object for animation
    states = mvs.nodes.States(snapshots=[
        mvs.nodes.Snapshot(state=state1),
        mvs.nodes.Snapshot(state=state2),
    ])

    # Save the multi-state object to an MVSJ file
    with open("5vpy_animation.mvsj", "w") as f:
        f.write(states.dumps())

    print("Saved animation MVSJ file to '5vpy_animation.mvsj'")


if __name__ == "__main__":
    create_animation()
