FROM heroku/miniconda

# Grab requirements.txt.
ADD ./Bioinformatic-App/requirements.txt /tmp/requirements.txt

# Install dependencies
RUN pip install -qr /tmp/requirements.txt

# Add our code
ADD ./Bioinformatic_App /opt/Bioinformatic_App/
WORKDIR /opt/Bioinformatic_App

RUN conda create -c rdkit -n my-rdkit-env rdkit

CMD gunicorn --bind 0.0.0.0:$PORT wsgi