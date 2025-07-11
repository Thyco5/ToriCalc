from rdkit import Chem
mol = Chem.MolFromSmiles("CCO") # Ethanol C0-C1-O2

# Test various CH2 SMARTS
smarts_ch2_v1 = "[CH2X4D2]"
smarts_ch2_v2 = "[CH2X4D2;!R]"
smarts_ch2_v3 = "[CH2X4D2;!R;!$(*=[O,S])]"
# This is the one that was puzzling, ensuring no other heteroatoms EXCEPT implicit oxygen (which is not in the list)
smarts_ch2_v4 = "[CH2X4D2;!R;!$(*~[#7,#9,#16,#17,#35,#53,F,Cl,Br,I]);!$(*=[O,S])]"
# Explicitly allow bonding to #6 (carbon) or #8 (oxygen) for its two non-H bonds,
# while still having the other general exclusions.
smarts_ch2_v5 = "[CH2X4D2;!R;!$(*=[O,S]);$([CH2X4D2]([#6,#8])[#6,#8]);!$(*~[#7,#9,#16,#17,#35,#53,F,Cl,Br,I])]"
# Let's simplify v5: CH2, not in ring, not carbonyl, and its neighbors are NOT (N, S, P, halogens)
# This means neighbors CAN BE C or O.
smarts_ch2_v6 = "[CH2X4D2;!R;!$(*=[O,S]);!$(*~[#7,#9,#16,#17,#35,#53,F,Cl,Br,I])]"


p_ch2_v1 = Chem.MolFromSmarts(smarts_ch2_v1)
p_ch2_v2 = Chem.MolFromSmarts(smarts_ch2_v2)
p_ch2_v3 = Chem.MolFromSmarts(smarts_ch2_v3)
p_ch2_v4 = Chem.MolFromSmarts(smarts_ch2_v4) # Your last attempt
p_ch2_v5 = Chem.MolFromSmarts(smarts_ch2_v5)
p_ch2_v6 = Chem.MolFromSmarts(smarts_ch2_v6)


print(f"Ethanol (CCO) C1 is CH2. Atom 0 (C), Atom 2 (O)")
print(f"CH2_SMARTS_V1 ({smarts_ch2_v1}): {mol.GetSubstructMatches(p_ch2_v1)}")
print(f"CH2_SMARTS_V2 ({smarts_ch2_v2}): {mol.GetSubstructMatches(p_ch2_v2)}")
print(f"CH2_SMARTS_V3 ({smarts_ch2_v3}): {mol.GetSubstructMatches(p_ch2_v3)}")
print(f"CH2_SMARTS_V4 ({smarts_ch2_v4}): {mol.GetSubstructMatches(p_ch2_v4)}") # This one gave ()
print(f"CH2_SMARTS_V5 ({smarts_ch2_v5}): {mol.GetSubstructMatches(p_ch2_v5)}")
print(f"CH2_SMARTS_V6 ({smarts_ch2_v6}): {mol.GetSubstructMatches(p_ch2_v6)}") # This is the same as v4

# For completeness, let's also make sure CH3 and OH still match as before
smarts_ch3 = "[CH3X4D1;!R;!$(*~[#7,#8,#9,#16,#17,#35,#53,F,Cl,Br,I]);!$(*=[O,S])]"
smarts_oh = "[OH1X2D1][#6X4;!c;!$(C=[O,N,S])]"
p_ch3 = Chem.MolFromSmarts(smarts_ch3)
p_oh = Chem.MolFromSmarts(smarts_oh)
print(f"CH3 SMARTS ({smarts_ch3}): {mol.GetSubstructMatches(p_ch3)}")
print(f"OH SMARTS ({smarts_oh}): {mol.GetSubstructMatches(p_oh)}")