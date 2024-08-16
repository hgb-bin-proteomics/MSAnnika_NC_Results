from fn import *

proteins = get_proteins()
interactions = get_interactions(proteins)
xi = get_xi()
print("---- xiSearch ----")
test_df(xi, interactions)
annika1 = get_annika_1()
print("---- MS Annika 1 ----")
test_df(annika1, interactions)
annika2 = get_annika_2()
print("---- MS Annika 2 ----")
test_df(annika2, interactions)
