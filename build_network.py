from bmtk.builder import NetworkBuilder
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list
import synapses

np.random.seed(123412)

# Initialize our network
net = NetworkBuilder("SPWR_biophysical")

# Create the possible x,y,z coordinates
xside_length = 600; yside_length = 600; height = 600; min_dist = 20;
x_grid = np.arange(0,xside_length+min_dist,min_dist)
y_grid = np.arange(0,yside_length+min_dist,min_dist)
z_grid = np.arange(0,height+min_dist,min_dist)
xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)
pos_list = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

#Number of cells in each population
numPN_A = 10800
numPN_C = 10800
numBask = 5400

# Load synapse dictionaries
synapses.load()
syn = synapses.syn_params_dicts()
###################################################################################
####################################Pyr Type A#####################################

# Pick coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numPN_A,replace=False)
pos = pos_list[inds,:]

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_A, pop_name='PyrA',
              positions=positions_list(positions=pos),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:feng_typeC',
              morphology=None)

##################################################################################
###################################Pyr Type C#####################################

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numPN_C,replace=False)
pos = pos_list[inds,:]

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_C, pop_name='PyrC',
              positions=positions_list(positions=pos),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:feng_typeC',
              morphology=None)
#################################################################################
############################# Chandelier ########################################

## Get rid of coordinates already used
#pos_list = np.delete(pos_list,inds,0)
#
## Pick new coordinates
#inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numAAC,replace=False)
#pos = pos_list[inds,:]
#
## Add a population of numAAC nodes
#net.add_nodes(N=numAAC, pop_name='AAC',
#              positions=positions_list(positions=pos),
#              mem_potential='e',
#              model_type='biophysical',
#              model_template='hoc:chandelier',
#              morphology=None)

#################################################################################
###########################Fast - spiking PV ints################################

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numBask,replace=False)
pos = pos_list[inds,:]

# Add a population of numBask nodes
net.add_nodes(N=numBask, pop_name='Bask',
              positions=positions_list(positions=pos),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:basket',
              morphology=None)
################################################################################
############################# BACKGROUND INPUTS ################################

# External inputs
thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(N=numPN_A+numPN_C,
                   pop_name='tON',
                   potential='exc',
                   model_type='virtual')

# External inputs
exc_bg_bask = NetworkBuilder('exc_bg_bask')
exc_bg_bask.add_nodes(N=numBask,
                   pop_name='tON',
                   potential='exc',
                   model_type='virtual')
##############################################################################
############################## CONNECT CELLS #################################

def dist_conn_perc(src, trg, min_dist=0.0, max_dist=300.0, min_syns=1, max_syns=2, A=0.2, B=0.2):
    
    sid = src.node_id
    tid = trg.node_id
    # No autapses
    if sid==tid:
        return None
    else:
        src_pos = src['positions']
        trg_pos = trg['positions']
    dist =np.sqrt((src_pos[0]-trg_pos[0])**2+(src_pos[1]-trg_pos[1])**2+(src_pos[2]-trg_pos[2])**2)
        #print("src_pos: {} trg_pos: {} dist: {}".format(src_pos,trg_pos,dist))        
    prob = A*np.exp(-B*dist)

    if dist <= max_dist and np.random.uniform() < prob:
        tmp_nsyn = np.random.randint(min_syns, max_syns)
        #print("creating {} synapse(s) between cell {} and {}".format(tmp_nsyn,sid,tid))
    else:
        tmp_nsyn = 0

    return tmp_nsyn

def one_to_one(source, target):
    
    sid = source.node_id
    tid = target.node_id
    if sid == tid:
    #print("connecting cell {} to {}".format(sid,tid))
        tmp_nsyn = 1
    else:
        return None

    return tmp_nsyn


def syn_dist_delay(source, target, min_delay):#, min_weight, max_weight):

	x_ind,y_ind,z_ind = 0,1,2

	dx = target['positions'][x_ind] - source['positions'][x_ind]
	dy = target['positions'][y_ind] - source['positions'][y_ind]
	dz = target['positions'][z_ind] - source['positions'][z_ind]

	dist = np.sqrt(dx**2 + dy**2 + dz**2)
	distDelay = dist/500
	#print("delay = {}".format(distDelay)) 
	return float(min_delay) + distDelay

def syn_dist_delay_section(source, target, min_delay, sec_id=0, sec_x=0.9):
    return syn_dist_delay(source, target, min_delay), sec_id, sec_x

# Create connections between Pyr --> Pyr cells
dynamics_file = 'PN2PN.json'

conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_one',
              connection_rule=dist_conn_perc,
              connection_params={'min_dist':0.0,'max_dist':300.0,
			         'min_syns':1,'max_syns':2,'A':0.01366,'B':0.008618},
              syn_weight=1,
	      delay=0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[0.0, 300.0],
              target_sections=['basal'])

#conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                    rule=syn_dist_delay_section,
#                    rule_params={'min_delay':syn[dynamics_file]['delay'], 
#				 'sec_id':0, 'sec_x':0.9},
#                    dtypes=[np.float, np.int32, np.float])

# Create connections between Pyr --> Bask cells
dynamics_file = 'PN2INT.json'

conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'Bask'},
              iterator = 'one_to_one',
	      connection_rule=dist_conn_perc,
              connection_params={'min_dist':0.0,'max_dist':300.0,
			         'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
              syn_weight=1,
	      delay = 0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[0.0, 300.0],
              target_sections=['somatic'])

#conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                    rule=syn_dist_delay_section,
#                    rule_params={'min_delay':syn[dynamics_file]['delay'], 
#				 'sec_id':0, 'sec_x':0.9},
#                    dtypes=[np.float, np.int32, np.float])



# Create connections between Pyr --> AAC cells
#net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'AAC'},
#              connection_rule=dist_conn_perc,
#              connection_params={'min_dist':0.0,'max_dist':300.0,
#			         'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
#              syn_weight=5.0e-03,
#              weight_function='lognormal',
#              weight_sigma=1.0e-03,
#              dynamics_params='AMPA_ExcToExc.json',
#              model_template='Exp2Syn',
#              distance_range=[0.0, 300.0],
#              target_sections=['somatic'],
#              delay=2.0)

# Create connections between Bask --> Pyr cells
dynamics_file = 'INT2PN.json'

conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_one',
              connection_rule=dist_conn_perc,
              connection_params={'min_dist':0.0,'max_dist':300.0,
			     'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
              syn_weight=1,
              delay=0.1,
              dynamics_params='INT2PN.json',
              model_template=syn['INT2PN.json']['level_of_detail'],
              distance_range=[0.0, 300.0],
              target_sections=['somatic'])

#conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                    rule=syn_dist_delay_section,
#                    rule_params={'min_delay':syn[dynamics_file]['delay'], 
#				 'sec_id':0, 'sec_x':0.9},
#                    dtypes=[np.float, np.int32, np.float])




# Create connections between AAC --> Pyr cells
#net.add_edges(source={'pop_name': 'AAC'}, target={'pop_name': ['PyrA','PyrC']},
#              connection_rule=dist_conn_perc,
#              connection_params={'min_dist':0.0,'max_dist':300.0,
#			     'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
#              syn_weight=5.0e-03,
#              weight_function='lognormal',
#              weight_sigma=1.0e-03,
#              dynamics_params='GABA_AAC.json',
#              model_template='Exp2Syn',
#              distance_range=[0.0, 300.0],
#              target_sections=['somatic'],
#              delay=2.0)

# Create connections between Bask --> Bask cells
dynamics_file = 'INT2INT.json'

conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['Bask']},
              iterator = 'one_to_one',
              connection_rule=dist_conn_perc,
              connection_params={'min_dist':0.0,'max_dist':300.0,
			     'min_syns':1,'max_syns':2,'A':0.24,'B':0.0},
              syn_weight=5.0e-03,
              weight_function='lognormal',
              weight_sigma=1.0e-03,
	      delay=0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[0.0, 300.0],
              target_sections=['somatic'])

#conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                    rule=syn_dist_delay_section,
#                    rule_params={'min_delay':syn[dynamics_file]['delay'], 
#				 'sec_id':0, 'sec_x':0.9},
#                    dtypes=[np.float, np.int32, np.float])

net.add_gap_junctions(source={'pop_name': ['Bask']}, 
		      target={'pop_name': ['Bask']},
 		      resistance = 0.0001, target_sections=['somatic'], 
		      connection_rule=dist_conn_perc,
		      connection_params={'min_dist':0.0,
					'max_dist':300.0,'min_syns':1,
					'max_syns':2,'A':0.08,'B':0.0})

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrA'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrC'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

net.add_edges(source=exc_bg_bask.nodes(), target=net.nodes(pop_name='Bask'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

print("Internal nodes and edges built")

# Create connections between "thalamus" and Pyramidals
# First define the connection rule

# Build and save our network

thalamus.build()
thalamus.save_nodes(output_dir='network')

exc_bg_bask.build()
exc_bg_bask.save_nodes(output_dir='network')
#
#print("External nodes and edges built")
t_sim = 5000

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=t_sim, dt = 0.1,
		report_vars = ['v'],
		spikes_inputs=[('mthalamus',   # Name of population which spikes will be generated for
                                'mthalamus_spikes.h5'),('exc_bg_bask','exc_bg_bask_spikes.h5')],
		#current_clamp={     
                #     'gids': [0],
                #     'amp': [0.5], 
                #     'delay': 100.0, 
                #     'duration': 50.0 
                # },
		components_dir='biophys_components',
		compile_mechanisms=True)


from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
#
psg = PoissonSpikeGenerator(population='mthalamus')
psg.add(node_ids=range(numPN_A+numPN_C),  # Have nodes to match mthalamus
        firing_rate=0.002,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim))    # Firing starts at 0 s up to 3 s
psg.to_sonata('mthalamus_spikes.h5')

psg = PoissonSpikeGenerator(population='exc_bg_bask')
psg.add(node_ids=range(numBask),  # Have nodes to match mthalamus
        firing_rate=0.004,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim))    # Firing starts at 0 s up to 3 s
psg.to_sonata('exc_bg_bask_spikes.h5')
