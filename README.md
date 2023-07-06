# OPED

OPED is a Optimized Prime Editing Design based on deep learning and transfer learning.

## Platform:

Linux 5.4.0-42-generic #46\~18.04.1-Ubuntu x86_64 GNU/Linux

## Conda env:

```bash
conda create -n OPED
conda activate OPED
conda install python=3.8.8 django=3.2.5  pytorch=1.6.0 cpuonly  pip sqlite pandas pymysql

conda install -c bioconda perl-dbi
conda install -c bioconda perl-dbd-sqlite
conda install -c anaconda virtualenv

conda install uwsgi

pip install biopython
```

## How to Run:

```bash
cd /your/path/to/OPED
conda activate OPED

# make soft link of OPED_nginx.conf into /etc/nginx/sites-enabled/
sudo ln -s OPED_nginx.conf /etc/nginx/sites-enabled/

# edit your nginx config file and copy it to your nginx config file location
# sudo cp /etc/nginx/nginx.conf /etc/nginx/nginx.conf.bakbeforeoped
sudo cp nginx.confg /etc/nginx/nginx.conf

# reload your nginx
nginx -s reload

# make migrations
python manage.py makemigrations
python manage.py migrate

# start the project
uwsgi --ini uwsgi.ini

# reload the project
uwsgi --reload OPED/OPED-master.pid

# stop the project
uwsgi --stop OPED/OPED-master.pid
```
