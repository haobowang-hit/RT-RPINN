from abaqus import session
from abaqusConstants import *
import csv
import os
import numpy as np
from collections import defaultdict

odb = session.odbs[session.odbs.keys()[-1]]
instance_name = odb.rootAssembly.instances.keys()[1]
instance = odb.rootAssembly.instances[instance_name]

for step_name, step in odb.steps.items():
    for frame_index, frame in enumerate(step.frames):

        displacement_field = frame.fieldOutputs['U']
        stress_field = frame.fieldOutputs['S']
        strain_field = frame.fieldOutputs['LE']

        displacement_values = displacement_field.getSubset(position=NODAL).values
        node_disp_map = {val.nodeLabel: val.data for val in displacement_values}

        stress_values = stress_field.getSubset(position=INTEGRATION_POINT).values
        node_stress_map = defaultdict(list)
        for value in stress_values:
            element = value.instance.elements[value.elementLabel-1]
            for node_label in element.connectivity:
                node_stress_map[node_label].append(value.data)
        average_node_stresses = {
            node_label: np.mean(stresses, axis=0) for node_label, stresses in node_stress_map.items()
        }

        strain_values = strain_field.getSubset(position=INTEGRATION_POINT).values
        node_strain_map = defaultdict(list)
        for value in strain_values:
            element = value.instance.elements[value.elementLabel-1]
            for node_label in element.connectivity:
                node_strain_map[node_label].append(value.data)
        average_node_strains = {
            node_label: np.mean(strains, axis=0) for node_label, strains in node_strain_map.items()
        }

        odb_basename = os.path.basename(odb.path).replace('.odb', '')
        output_filename = '{}_{}_frame{}.csv'.format(odb_basename, step_name, frame_index)

        with open(output_filename, mode='wb') as file:
            writer = csv.writer(file)

            writer.writerow([
                'NodeLabel', 'X', 'Y', 'Z', 'U1', 'U2', 'U3',
                'LE11', 'LE22', 'LE33', 'LE12', 'LE13', 'LE23',
                'S11', 'S22', 'S33', 'S12', 'S13', 'S23'
            ])

            for node in instance.nodes:
                node_label = node.label
                X, Y, Z = node.coordinates

                U = node_disp_map.get(node_label, [0.0, 0.0, 0.0])
                strain_data = average_node_strains.get(node_label, [0.0]*6)
                LE11, LE22, LE33, LE12, LE13, LE23 = strain_data
                stress_data = average_node_stresses.get(node_label, [0.0]*6)
                S11, S22, S33, S12, S13, S23 = stress_data

                writer.writerow([
                    node_label, X, Y, Z, U[0], U[1], U[2],
                    LE11, LE22, LE33, LE12, LE13, LE23,
                    S11, S22, S33, S12, S13, S23
                ])

        print('Saved: {}'.format(output_filename))
