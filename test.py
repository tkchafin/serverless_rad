import os 
import sys 

from lithopsrad.manager import PipelineManager

# NOTE: This pipeline is currently intended as a demonstration, 
# and thus does not currently have advanced features such as checkpointing 
# or orchestration 

def main():

    # get params for lithops and run
    config_file = "lithops_config"
    
    # configure pipeline 
    workflow = PipelineManager(config_file)

    # run pipeline 
    workflow.run()

if __name__ == "__main__":
    main()