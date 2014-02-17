INSTALL

virtualenv venv
source venv/bin/activate
pip install -r dependencies.txt

RUN

virtualenv venv
ipython notebook --pylab inline
