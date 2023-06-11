import pandas as pd
import numpy as np


class Data(object):
    def __init__(self, virus='HCV', vaccine_status=False, sheet_name='sheet1'):
        self.filename = './VL_data.xlsx'
        self.virus = virus
        self.vaccine_status = vaccine_status
        self.sheet_name = sheet_name
        self.load_data()

    def load_data(self):
        df = pd.read_excel(self.filename, sheet_name=self.sheet_name)

        # select data
        data = df[df['Virus'] == self.virus]
        data = data[data['Vaccinated'] == self.vaccine_status]
        patient_seq = np.unique(np.array(data['Patient']))

        output = {'id': [], 'lgVL': [], 'time': [],
                  'symptom': [], 'LOD': []}
        for num_patient in range(len(patient_seq)):
            # --------Fit Viral Loads--------------------------
            data_tmp = data[data['Patient'] == patient_seq[num_patient]]
            data_tmp = data_tmp.dropna(axis=0, how='all', subset=['lgVL', ])  # drop lgVL=nan
            t_point = np.array(data_tmp['Time'], dtype=np.int32)
            lgVL = np.array(data_tmp['lgVL'])

            output['id'].append(patient_seq[num_patient])
            output['lgVL'].append(lgVL)
            output['time'].append(t_point)
            output['symptom'].append(np.array(data_tmp['Symptom'])[0])
            output['LOD'].append(data_tmp['LOD'].iloc[0])

        self.data = output

    def single_patient_data(self, id):
        if id < len(self.data['id']):
            return {'id': self.data['id'][id],
                    'lgVL': self.data['lgVL'][id],
                    'time': self.data['time'][id],
                    'symptom': self.data['symptom'][id],
                    'LOD': self.data['LOD'][id]}
        else:
            print('Please select between [%s]!' % range(len(self.data['id'])))
