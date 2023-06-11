import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class Data(object):
    def __init__(self):
        
        self.filename = './data_Tan2021.xlsx'
        self.load_data()

    def load_data(self):
        
        df_VL = pd.read_excel(self.filename, sheet_name='VL')

        df_Imm = pd.read_excel(self.filename, sheet_name='Imm')

        id_seq = np.unique(np.array(df_VL['ID']))

        output = {'id': [], 'group': [],
                  'lgVL': [], 't_VL': [], 'lgLOD': [],
                  'T cell': [], 'Ab': [], 't_Imm': [],
                  'Ab x T': []
                  }
        
        for num in range(len(id_seq)):
            
            # VL
            data_tmp = df_VL[df_VL['ID'] == id_seq[num]]

            t_point = np.array(data_tmp['Days'], dtype=np.int32)
            lgVL = np.array(data_tmp['lgVL'])

            output['id'].append(id_seq[num])
            output['group'].append(data_tmp['Group'].iloc[0])
            output['lgVL'].append(lgVL)
            output['t_VL'].append(t_point)
            output['lgLOD'].append(data_tmp['lgLOD'].iloc[0])

            # T cell and Ab
            data_tmp = df_Imm[df_Imm['ID'] == id_seq[num]]
            output['T cell'].append(np.array(data_tmp['T cell']))
            output['Ab'].append(np.array(data_tmp['Ab']))
            output['Ab x T'].append(np.array(data_tmp['Ab'])*np.array(data_tmp['T cell']))
            output['t_Imm'].append(np.array(data_tmp['Days'], dtype=np.int32))


        self.dataset = output


    def single_patient_data(self, id):
        if id in self.dataset['id']:
            num = self.dataset['id'].index(id)
            return {k:v[num]
                    for k, v in zip(self.dataset.keys(), self.dataset.values())}
        else:
            print('Please select in following id: \n', self.dataset['id'])

    
    def plot_time_course(self, t_by_peak=False):
        
        fig, ax = plt.subplots(2,2,figsize=(8,7))
        ax=ax.flat

        for id in self.dataset['id']:
            dt = self.single_patient_data(id)

            if dt['group'] == 'Severe':
                clr='r'
            elif dt['group'] == 'Moderate':
                clr='orange'
            else:
                clr='green'

            if t_by_peak:
                peak_t = dt['t_VL'][np.nanargmax(dt['lgVL'])]
                dt['t_VL'] = dt['t_VL'] - peak_t
                dt['t_Imm'] = dt['t_Imm'] - peak_t


            ax[0].plot(dt['t_VL'], dt['lgVL'], 
                       c=clr, lw=0.8, alpha=0.5,
                       marker='o', markersize=4)
            ax[1].plot(dt['t_Imm'], dt['T cell'], c=clr, lw=0.8)
            ax[2].plot(dt['t_Imm'], dt['Ab'], c=clr, lw=0.8)
            ax[3].plot(dt['t_Imm'], dt['Ab']*dt['T cell'], c=clr, lw=0.8)

        if t_by_peak:
            xlab = 'Days post peak VL'
            xlim = [-20, 20]
        else:
            xlab = 'Days post symptom onset'
            xlim = [0, 50]

        for a, ylab in zip(ax, ['lgVL', 'T cell', 'Ab inhibition (%)', 'T x Ab']):
            a.set_ylabel(ylab)
            a.set_xlim(xlim)
            a.set_xlabel(xlab)
        
        ax[0].set_ylim([0.5, 8])
        ax[1].set_ylim([-10, 3200])
        ax[2].set_ylim([-5, 100])
        fig.tight_layout()




if __name__ == '__main__':
    d = Data()
    d.plot_time_course(t_by_peak=True)
