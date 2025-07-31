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
    - title 'State 2'
    - description: 'The switch II region' 
    - highlight the residues 22 - 44 in a different colour. 
    """


import molviewspec as mvs

# 1. Create a builder
builder = mvs.create_builder()

# 2. Download and parse the structure from PDB (mmCIF format)
structure = builder.download(url="https://files.rcsb.org/download/5vpy.cif").parse(format="mmcif").model_structure()

# 3. Add representations for different components of the structure
# Show protein as a surface
structure.component(selector="protein").representation(type="surface")

# Show ligands as ball-and-stick
structure.component(selector="ligand").representation(type="ball_and_stick")

# 4. Save the state to an MVSJ file
builder.save_state(destination="5vpy.mvsj")

print("Saved MVSJ file to '5vpy.mvsj'")
