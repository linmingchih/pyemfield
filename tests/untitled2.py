from pyaedt import Hfss
import os

def export_ffd(output_folder, solution='', freq='', version='2024.2'):
    hfss = Hfss(version=version)
    ports = hfss.get_all_sources()
    
    info = {}
    for sol in hfss.post.available_report_solutions("Far Fields"):
        data = hfss.post.get_solution_data_per_variation('Far Field', setup_sweep_name=solution,  expressions='')
        unit = data.units_sweeps['Freq']
        freqs = [f'{f}{unit}' for f in data.primary_sweep_values]    
        info[sol] = freqs
    
    assert info, "No any solution exist!"    
    print('Existed solution:', info)
    
    if not solution:
        solution, freqs = info[info.keys()[0]]
        freq = freqs[0]
    
    if freq in info[solution]:
        print(f'{solution}{freq}')
    else:
        raise Exception(f'{solution}{freq} does not exist!')
    
    x = hfss.insert_infinite_sphere(x_start=0, x_stop=360, x_step=1, y_start=0, y_stop=180, y_step=1)
    data = hfss.post.get_solution_data_per_variation('Far Field', setup_sweep_name=solution,  expressions='')
    unit = data.units_sweeps['Freq']
    freqs = [f'{f}{unit}' for f in data.primary_sweep_values]    
    
    if not freq:
        freq = f'{freqs[0]}'
    
    oModule = hfss.odesign.GetModule("RadField")
    try:
        for p1 in ports:
            setting = {}
            for p2 in ports:
                if p2 == p1:
                    setting[p2] = ('1W', '0deg')
                else:
                    setting[p2] = ('0W', '0deg')
            
            hfss.edit_sources(setting)
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

export_ffd('d:/demo2', 'Setup1 : Sweep', '27.5GHz')
