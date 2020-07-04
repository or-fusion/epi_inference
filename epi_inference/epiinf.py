import os
from . import workflow
from epi_inference.engine.driver import driver
  
def main():
    if "HOME" in os.environ:
        tmpdir = os.path.join(os.environ["HOME"],".epiinf","tmp")
    else:
        tmpdir = os.path.join(os.getcwd(),"epiinf_tmp")
    driver(tmpdir=tmpdir, command='epiinf')

if __name__ == "__main__":
    main()
