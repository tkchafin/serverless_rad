import os 
import sys 

from lithopsrad.manager import PipelineManager

def main():

    # get params for lithops and run
    config_file = "lithops_config"
    
    workflow = PipelineManager(config_file)
    workflow.run()

if __name__ == "__main__":
    main()