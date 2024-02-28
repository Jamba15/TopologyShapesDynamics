from functions_signalprocess import *
config = config_parser('signal_processing_config.yaml')
graph = SimplicialComplexConstructor(**config['simplicial_complex_constructor'])
graph.draw_graph()
#%%
# Signal Processing on dim 0-1, using D1
if config['process_with'] == 'D1':
    _, D = graph.create_B1_D1()
else:
    _, D = graph.create_B2_D2()

signals = signal_elaboration(D)
signals.plot_results(save_fig=False, fig_name='Processing with '+config['process_with'])

