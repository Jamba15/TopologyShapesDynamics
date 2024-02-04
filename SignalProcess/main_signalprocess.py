from functions_signalprocess import *

multiplex_name = "Padgett-Florentine-Families_multiplex.edges"
edge_list = "NGF_edgelist.txt"
trinangles_list = "NGF_triangles_d2.txt"

graph = GraphCreator(edge_list_file=edge_list, triangle_file=trinangles_list)
# graph.draw_graph()
#%%
B1, D1 = graph.create_B1_D1()
B2, D2 = graph.create_B2_D2()

# Signal Processing on dim 0-1, using D1
D = D1.copy()

signals = signal_elaboration(D)
signals.plot_results(save_fig=True, fig_name="D1_signals")
