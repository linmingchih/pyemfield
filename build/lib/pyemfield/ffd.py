import os, re
import numpy as np
from cmath import exp, sqrt
from scipy.interpolate import interp1d
from itertools import product
from math import sin, cos, radians, pi, log10
import matplotlib.pyplot as plt
from collections import defaultdict
from pyaedt import Hfss
from scipy.optimize import minimize



def dbm_to_w(x):
    return pow(10, x / 10) / 1000


def calculate_cell_size():
    data = []
    for theta, phi in product(range(181), range(361)):
        data.append(sin(radians(theta))*radians(1)*radians(1))
    
    return np.array(data).reshape(181, 361) / sum(data) * 4 * pi

class Cdf:
    def __init__(self, name, x, y, xlabel, title=None):        
        self.data = [(name, x, y)]
        plt.title(name)
        plt.xlabel(xlabel)
        plt.ylabel('Cumulative Probability')
        if title:
            plt.title(title)

    def __add__(self, other):
        self.data += other.data
        return self
    
    def plot(self):
        colors = ['blue','green','red','cyan','magenta','yellow','black']
        for (name, x, y), c in zip(self.data, colors):
            plt.plot(x, y, label=name, color=c)
        plt.grid()
        plt.legend()
        plt.show()


class Ffd:
    cell_size = calculate_cell_size()
    
    def __init__(self, ffd_path='', ffd_name=''):
        self.name = ffd_name
        
        self._etheta = np.zeros((181, 361))
        self._ephi = np.zeros((181, 361))
        self.radiated_power_density = np.zeros((181, 361))
        
        if ffd_path:
            assert os.path.isfile(ffd_path), f"{ffd_path} is not valid."        
            self.etheta, self.ephi = self.get_ffd(ffd_path)        
        
            print(ffd_name)
        

    @property
    def etheta(self):
        return self._etheta
    
    @etheta.setter
    def etheta(self, value):
        self._etheta = value
        self.radiated_power_density = (np.power(np.absolute(self.etheta),2)
                               + np.power(np.absolute(self.ephi),2))/377/2
        self._sum_radiation = np.sum(self.radiated_power_density * Ffd.cell_size)
        
        
    @property
    def ephi(self):
        return self._ephi
    
    @ephi.setter
    def ephi(self, value):
        self._ephi = value
        self.radiated_power_density = (np.power(np.absolute(self.etheta),2)
                               + np.power(np.absolute(self.ephi),2))/377/2
        self._sum_radiation = np.sum(self.radiated_power_density * self.cell_size)
        
    @staticmethod
    def get_ffd(ffd_path):
        with open(ffd_path) as f:
            text = f.readlines()
        
        assert text[0] == '0 180 181\n', "invalid ffd theta resolution!"
        assert text[1] == '0 360 361\n', "invalid ffd phi resolution!"
        
        etheta = []
        ephi = []
        
        for line in text[4:]:
            etheta_re, etheta_im, ephi_re, ephi_im = map(float, line.strip().split())
            etheta.append(complex(etheta_re, etheta_im))
            ephi.append(complex(ephi_re, ephi_im))
        
        return np.reshape(np.array(etheta), (181, 361)), np.reshape(np.array(ephi), (181, 361))        
    
    @staticmethod
    def get_loss(s2p_path):
        mapping = {'hz':1, 'khz':1e3, 'mhz':1e6, 'ghz':1e9, 'thz':1e12}
        
        with open(s2p_path) as f:
            text = f.readlines()
        
        data = []
        for line in text:
            try:
                data += [float(i) for i in line.strip().split()]
                continue
            except:
                pass
            
            m = re.search('#\s+(\w+)\s+S\s+MA\s+R\s+50', line)
            if m:
                scale = mapping[m.group(1).lower()]
            
        
        n = 0
        for freq in data[::9]:
            if freq*scale == Ffd.frequency:
                return (data[n+3], data[n+4])
            n+=9

    def __add__(self, other):
        result = Ffd()
        result.etheta = self.etheta + other.etheta
        result.ephi = self.ephi + other.ephi

        return result
    
    def set(self, mag, phase):
        result = Ffd()
        result.etheta = np.sqrt(mag)*np.exp(1j*np.radians(phase))*self.etheta
        result.ephi = np.sqrt(mag)*np.exp(1j*np.radians(phase))*self.ephi
        
        return result    
            
    def plot_rETotal(self):
        self.retotal = np.sqrt(np.power(np.absolute(self.etheta),2)+np.power(np.absolute(self.ephi),2))
  
        plt.figure(figsize=(8, 4))
        plt.title(f'{self.name} - rETotal')
        plt.xlabel('Phi (deg)')
        plt.ylabel('Theta (deg)')
        
        plt.imshow(self.retotal, cmap='jet')
        
        plt.colorbar()
        
    def __repr__(self):
        return self.name

class Beam:
    def __init__(self, *argc):
        if len(argc)==1 and type(argc[0]==dict):
            ffd_excitation = argc[0]
            self.ffds = ffd_excitation.keys()
            self.ffd_excitation = ffd_excitation
            self.set_excitations(self.ffd_excitation)
            
        elif len(argc)==3 and type(argc[0])==list:
            self.ffds, mag, phase = argc[0], argc[1], argc[2]
            self.ffd_excitation = {i:(mag, phase) for i in self.ffds}
            self.set_excitations(self.ffd_excitation)
        
        else:
            raise Exception('Error: Beam initialization is invalid!')
    
    # def __repr__(self):
    #     return '\n'.join([f'{i.name}:({mag}, {phase})' for i, (mag, phase) in self.ffd_excitation.items()])
    
    def set_excitations(self, ffd_excitation):
        self.ffd_excitation = ffd_excitation
        self._sum_excitation = 0
        self._sum_ffd = Ffd()
        
        for i, (mag, phase) in self.ffd_excitation.items():
            self._sum_ffd += i.set(mag, phase)
            self._sum_excitation += mag
        
        
        self._sum_radiation = np.sum(self._sum_ffd.radiated_power_density)
        self.realized_gain = self._sum_ffd.radiated_power_density/(self._sum_excitation/4/np.pi)
    
    def plot_realized_gain_contour(self, name=''):
        plt.figure(figsize=(8, 4))
        if name:
            plt.title(f'{name}: Realized Gain(dB)')
        else:
            plt.title('Realized Gain(dB)')
        plt.xlabel('Phi (deg)')
        plt.ylabel('Theta (deg)')
        
        plt.imshow(10*np.log10(self.realized_gain), cmap='jet')
        plt.colorbar()
        plt.contour(10*np.log10(self.realized_gain), colors='black', linewidths=0.5)
        plt.show()
        print(10*log10(np.max(self.realized_gain)))

    def _get_realized_gain_cdf(self, name=''):
        cell_size = calculate_cell_size()
        
        data = []
        for theta in range(181):
            for phi in range(361):
                data.append((self.realized_gain[theta, phi], cell_size[theta, phi]))

        
        data = sorted(data)
        x, y = [], []
        accumulation = 0
        for realized_gain, area in data:
            accumulation += area
            x.append(10*log10(realized_gain))
            y.append(accumulation/4/pi)

        f = interp1d(y, x)
        y_new = [1, 0.8, 0.5, 0.2]
        x_new = f(y_new)
        for u, v in zip(y_new, x_new):
            print(f'{u*100}% CDF is {v:.3f}')      
        
        return Cdf(name, x, y, 'Realized Gain')
        
    def plot_realized_gain_cdf(self, name=''):
        cdf = self._get_realized_gain_cdf(name='')
        cdf.plot()
        return cdf        
        
    
    def get_realized_gain(self, theta, phi):
        return self.realized_gain[theta, phi]
    

    def optimize_gain_exhaustive(self, theta, phi, bits=4):
        etheta = [i._etheta[theta, phi]*m for i, (m, p) in self.ffd_excitation.items()]
        ephi = [i._ephi[theta, phi]*m for i, (m, p) in self.ffd_excitation.items()]
        
        phase_list = [[0]]+[[(360/(2**bits))*j for j in range(2**bits)] for i in range(len(etheta)-1)]

        et_list = []
        for phase , et in zip(phase_list, etheta):
            et_list.append([et*exp(1j*radians(p)) for p in phase])
        
        ep_list = []    
        for phase , ep in zip(phase_list, ephi):
            ep_list.append([ep*exp(1j*radians(p)) for p in phase])
        
        max_mag = float('-inf')
        max_phase = None
        
        for p, et, ep in zip(product(*phase_list), product(*et_list), product(*ep_list)):
            mag = sqrt(abs(sum(et))**2 + abs(sum(ep))**2)
            if mag.real > max_mag:
                max_mag = mag.real
                max_phase = p

        result = {}
        for (i, (mag, phase)), j in zip(self.ffd_excitation.items(), max_phase):
            result[i] = (mag, j)
        
        b = Beam(result)
        b.phase_code = [int(i/(360/(2**bits))) for i in max_phase]
        return b

    def optimize_gain_heuristic(self, theta, phi, bits=4):
        etheta = [i._etheta[theta, phi]*m for i, (m, p) in self.ffd_excitation.items()]
        ephi = [i._ephi[theta, phi]*m for i, (m, p) in self.ffd_excitation.items()]
        phases = [(360/(2**bits))*i for i in range(2**bits)]
        
        et0 = etheta[0]
        ep0 = ephi[0]
        et_list = [[et0]]
        ep_list = [[ep0]]
        phase_list = [[0]]
        for et, ep in zip(etheta[1:], ephi[1:]):
            data = []
            for p in phases:
                mag = abs(et*exp(1j*radians(p)) + et0)**2 + abs(ep*exp(1j*radians(p)) + ep0)**2
                data.append((mag.real, et, ep, p))
            data = sorted(data, reverse=True)
            m, x_et, x_ep, x_phase = zip(*data[0:2**(bits-1)])
            et_list.append(x_et)
            ep_list.append(x_ep)
            phase_list.append(x_phase)

        max_mag = float('-inf')
        max_phase = None        
        for p, et, ep in zip(product(*phase_list), product(*et_list), product(*ep_list)):
            mag = sqrt(abs(sum(et))**2 + abs(sum(ep))**2)
            if mag.real > max_mag:
                max_mag = mag.real
                max_phase = p
                
        result = {}
        for (i, (mag, phase)), j in zip(self.ffd_excitation.items(), max_phase):
            result[i] = (mag, j)
        
        b = Beam(result)
        b.phase_code = [int(i/(360/(2**bits))) for i in max_phase]            
        return b



    def optimize_gain(self, theta, phi):
        # 提取 etheta 和 ephi 值
        etheta = [i._etheta[theta, phi] * m for i, (m, p) in self.ffd_excitation.items()]
        ephi = [i._ephi[theta, phi] * m for i, (m, p) in self.ffd_excitation.items()]

        # 定義目標函數（負的增益值，因為 minimize 是尋找最小值）
        def objective(phases):
            et_sum = sum(et * exp(1j * radians(p)) for et, p in zip(etheta, phases))
            ep_sum = sum(ep * exp(1j * radians(p)) for ep, p in zip(ephi, phases))
            gain_magnitude = abs(et_sum)**2 + abs(ep_sum)**2
            return -gain_magnitude.real  # 最大化增益 -> 最小化負增益

        # 初始位相設定為 0 度
        initial_phases = [0] * len(etheta)

        # 使用 scipy.optimize.minimize 進行優化
        result = minimize(objective, initial_phases, bounds=[(0, 360)] * len(etheta))

        # 獲取最佳位相
        optimized_phases = result.x

        # 將最佳位相應用於 Beam 物件
        result_dict = {}
        for (i, (mag, phase)), p in zip(self.ffd_excitation.items(), optimized_phases):
            result_dict[i] = (mag, p)

        # 更新 Beam 物件
        self.phase_code = optimized_phases  # 儲存最佳位相值
        return Beam(result_dict)  # 返回新的 Beam 物件，包含優化結果



class Plane:
    def __init__(self, name='', beams=[]):
        self.name = name
        self.beams = beams
        self.radiated_power_density = np.zeros((181, 361))
        for b in self.beams:
            self.radiated_power_density = np.maximum(self.radiated_power_density, b.sum_ffd.radiated_power_density)

        self.realized_gain = np.zeros((181, 361))
        for b in self.beams:
            self.realized_gain = np.maximum(self.realized_gain, b.realized_gain)  
        
        self.radiated_power = self.radiated_power_density*4*pi
        
    def __add__(self, other):
        name = self.name + '+' + other.name
        beams = self.beams + other.beams
        return Plane(name, beams)
    
    def get_rGain_cdf(self):
        print(f'\n-------------------{self.name} rGain------------------')
        cell_size = calculate_cell_size()
        data = []
        for theta in range(181):
            for phi in range(361):
                if self.realized_gain[theta, phi] > 0:
                    data.append((self.realized_gain[theta, phi], cell_size[theta, phi]))        

        data = sorted(data)
        x, y = [], []
        accumulation = 0
        for power, area in data:
            accumulation += area
            x.append(10*log10(power))
            y.append(accumulation/4/pi)
        
        f = interp1d(y, x)
        y_new = [1, 0.8, 0.5, 0.2]
        x_new = f(y_new)
        for u, v in zip(y_new, x_new):
            print(f'{u*100}% rGain CDF is {v:.3f}(dB)')        
        
        return Cdf(self.name, x, y, 'rGain')
    
    def plot_rGain_cdf(self):
        cdf = self.get_rGain_cdf()
        cdf.plot()
        return cdf
    
    def get_eirp_cdf(self):
        print(f'\n-------------------{self.name} EIRP-------------------')
        cell_size = calculate_cell_size()
        self.radiate_power = self.radiated_power_density*4*pi

        total_radiate_power = np.sum(self.radiated_power*cell_size)
        
        data = []
        for theta in range(181):
            for phi in range(361):
                if self.radiate_power[theta, phi] > 0:
                    data.append((self.radiated_power[theta, phi], cell_size[theta, phi]))

        
        data = sorted(data)
        x, y = [], []
        accumulation = 0
        for power, area in data:
            accumulation += area
            x.append(10*log10(power)+30)
            y.append(accumulation/4/pi)
        
        f = interp1d(y, x)
        y_new = [1, 0.8, 0.5, 0.2]
        x_new = f(y_new)
        for u, v in zip(y_new, x_new):
            print(f'{u*100}% EIRP CDF is {v:.3f}(dBm)')        
        
        return Cdf(self.name, x, y, 'EIRP (dbm)')
        
    def plot_eirp_cdf(self):
        self.get_eirp_cdf().plot()
    
    def plot_eirp_contour(self):
        plt.figure(figsize=(8, 4))
        plt.title('radiated power density(dBm)')
        plt.xlabel('Phi (deg)')
        plt.ylabel('Theta (deg)')
        
        plt.imshow(10*np.log10(self.radiated_power)+30, cmap='jet')
        plt.colorbar()
        plt.contour(10*np.log10(self.radiated_power)+30, colors='black', linewidths=0.5)
        
        plt.show()
        
    def plot_rGain_contour(self):
        plt.figure(figsize=(8, 4))
        plt.title('realized gain (dB)')
        plt.xlabel('Phi (deg)')
        plt.ylabel('Theta (deg)')
        
        plt.imshow(10*np.log10(self.realized_gain), cmap='jet')
        plt.colorbar()
        plt.contour(10*np.log10(self.realized_gain), colors='black', linewidths=0.5)
        
        plt.show()       
        
def create_plane(plane_name, _beam, angles, fast=True, bits=4):
    beams = []
    print(f'\n------------------Computing-{plane_name}-------------------')
    
    for theta, phi in angles:
        if fast:
            x = _beam.realized_gain_optimize2(theta, phi, bits)
        else:
            x = _beam.realized_gain_optimize(theta, phi, bits)
        beams.append(x)

        phases = [round(i[1],1) for i in x.ffd_excitation.values()]
        rgain = round(10*log10(x.get_realized_gain(theta, phi)), 2)
        peak = round(10*log10(np.max(x.realized_gain)), 2)
        print(f'To {theta},{phi}: Phase={phases}, rGain={rgain}(dB) Peak={peak}(dB)')
    return Plane(plane_name, beams)


def get_obj_name(x = locals()):
    result = {}
    for i, j in x.items():
        if type(j) in [Plane, Beam, Ffd]:
            result[j] = (i)
    return result

    
def dump_code(code_path, planes, mapping_table):
    
    module_num_dic = {'A':1, 'B':0}
    obj_name_dic = get_obj_name()
    
    text = ['\t'.join(['TRx', 'Pol', 'Band', 'Beam_ID', 'Module', 'Ant_Type', 'Beam_Type', 'NV', 'Phase'])]

    for plane in planes:
        beam_id = 0
        plane_name = obj_name_dic[plane]
        m = re.search('(A|B)plane_(HandV|HorV|H|V)', plane_name)
        if m:
            module_id = module_num_dic[m.group(1)]
            polarization = m.group(2)
        else:
            continue
        for beam in plane.beams:            
            ffd_id = f"({'-'.join([str(mapping_table.index(i)) for i in beam.ffds])})"
            phases = f"({'-'.join(map(str, beam.phase_code))})"
            line = ['Tx', polarization, Ffd.frequency/1e9, beam_id, module_id, 'Patch', 'NARROW', ffd_id, phases]
            text.append('\t'.join(map(str, line)))
            beam_id += 1
    
    with open(code_path, 'w') as f:
        f.write('\n'.join(text))

def load_code(code_path, mapping_table, mag=1, bits=4):
    group = defaultdict(list)
    result = {}
    
    with open(code_path) as f:
        text = f.readlines()

    for line in text[1:]:
        try:
            _, polaization, _, beam_id, module_id, _, _, ffd_ids, phases = line.strip().split('\t')
            module = 'A' if module_id == '1' else 'B'
            ffds = [mapping_table[int(i)] for i in ffd_ids.replace('(','').replace(')', '').split('-')]
            phases = [int(i)*(360/2**bits) for i in phases.replace('(','').replace(')', '').split('-')]
            beam = Beam({ffd:(mag, phase) for ffd, phase in zip(ffds, phases)})        
            group[(polaization, module)].append(beam)
        
        except:
            raise Exception(f'Error:{line}')

    for (polarization, module), beams in group.items():
        name = f'{module}plane_{polarization}'
        result[name] = Plane(name, beams)
    
    return result


class hfss_design:
    def __init__(self, version='2024.2'):
        self.hfss = Hfss(version=version)
        
        self.info = {}
        for sol in self.hfss.post.available_report_solutions("Far Fields"):
            data = self.hfss.post.get_solution_data_per_variation('Far Field', setup_sweep_name=sol,  expressions='')
            unit = data.units_sweeps['Freq']
            freqs = [f'{f}{unit}' for f in data.primary_sweep_values]    
            self.info[sol] = freqs
            
        self.ports = self.hfss.get_all_sources()
        print(f'Ports:{self.ports}')
        print('Existed solution:', self.info)    
        
    def update_excitation(self, excitations):
        excitations = {i.name:j for i, j in excitations.items()}
        sources = {}
        for port in self.ports:
            if port in excitations:
                x, y = excitations[port]
                sources[port] = (f'{x}W', f'{y}deg')
            else:
                sources[port] = ('1W', '0deg')
                
        self.hfss.edit_sources(sources)
        
        
    def export_ffds(self, output_folder, solution='', freq=''):
        assert self.info, "No any solution exist!"   

        if not solution:
            solution, freqs = self.info[self.info.keys()[0]]
            freq = freqs[0]
        
        if freq in self.info[solution]:
            print(f'{solution}{freq}')
        else:
            raise Exception(f'{solution}{freq} does not exist!')
        
        x = self.hfss.insert_infinite_sphere(x_start=0, x_stop=180, x_step=1, y_start=0, y_stop=360, y_step=1)
        data = self.hfss.post.get_solution_data_per_variation('Far Field', setup_sweep_name=solution,  expressions='')
        unit = data.units_sweeps['Freq']
        freqs = [f'{f}{unit}' for f in data.primary_sweep_values]    
        
        if not freq:
            freq = f'{freqs[0]}'
        
        oModule = self.hfss.odesign.GetModule("RadField")
        try:
            for p1 in self.ports:
                setting = {}
                for p2 in self.ports:
                    if p2 == p1:
                        setting[p2] = ('1W', '0deg')
                    else:
                        setting[p2] = ('0W', '0deg')
                
                self.hfss.edit_sources(setting)
                ffd_path = os.path.join(output_folder, f'{p1}.ffd')
            
                oModule.ExportFieldsToFile(
               	[
               		"ExportFileName:="	, ffd_path,
               		"SetupName:="		, x.name,
               		"IntrinsicVariationKey:=", f"Freq=\'{freq}\'",
               		"DesignVariationKey:="	, "",
               		"SolutionName:="	, solution,
               		"Quantity:="		, ""
               	])
                print(solution, freq, ffd_path)
        except:
            raise Exception('Export Error!')
        finally:
            oModule.DeleteSetup([x.name])
        
def get_ffds(folder):
    result = {}
    
    for ffd in os.listdir(folder):
        if not ffd.endswith('.ffd'):
            continue
        
        ffd_path = os.path.join(folder, ffd)
        key = ffd.replace('.ffd', '')
        ffd_obj = Ffd(ffd_path, key)
        result[key] = ffd_obj
        
    return result

if __name__ == '__main__':
    hd = hfss_design()
    
    folder = r'D:\OneDrive - ANSYS, Inc\GitHub\pyemf\tests\ffds'
    #hd.export_ffds(folder, 'Setup1 : Sweep', '30.0GHz')
    
    ffds = get_ffds(folder)
    
    ffds.keys()
    x = {j:(1,0) for i, j in ffds.items()}
    b1 = Beam(x)
    
    b1.ffd_excitation
    b2 = b1.optimize_gain(40, 60)
    b2.plot_realized_gain_contour()
    b2.ffd_excitation
    
