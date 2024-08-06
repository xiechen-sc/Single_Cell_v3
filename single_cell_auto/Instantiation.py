from scClass import *



def get_script_fun(module_analysis,config_path,yaml_data,project_id):
        if module_analysis == 'featureplot':
            sc_obj = Featureplot(project_id=project_id,config_path=config_path,**yaml_data)
        
        elif module_analysis == 'modified_cell_type':
            sc_obj = Modified_cell_type(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'diff':
            sc_obj = Diff(project_id=project_id,config_path=config_path,**yaml_data)
            
        sc_obj.get_script()