from scClass import *



def get_script_fun(module_analysis,config_path,yaml_data,project_id):
        if module_analysis == 'featureplot':
            sc_obj = Featureplot(project_id=project_id,config_path=config_path,**yaml_data)
        
        elif module_analysis == 'modified_cell_type':
            sc_obj = Modified_cell_type(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'diff':
            sc_obj = Diff(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'sub_clusters':
            sc_obj = Sub_Clusters(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'sub_clusters_old':
            sc_obj = Sub_Clusters_old(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'singleR':
            sc_obj = SingleR(project_id=project_id,config_path=config_path,**yaml_data)
        
        elif module_analysis == 'enrichment':
            sc_obj = Enrichment(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'scenic':
            sc_obj = Scenic(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'decontX':
            sc_obj = DecontX(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'monocle2':
            sc_obj = Monocle2(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'addmodulescore':
            sc_obj = Addmodulescore(project_id=project_id,config_path=config_path,**yaml_data)

        elif module_analysis == 'cellchat':
            sc_obj =  Cellchat(project_id=project_id,config_path=config_path,**yaml_data)
            
        sc_obj.get_script()
