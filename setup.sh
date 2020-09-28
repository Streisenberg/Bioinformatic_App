mkdir -p ~/.streamlit/
conda create -c rdkit -n my-rdkit-env rdkit
conda activate my-rdkit-env

echo "\
[server]\n\
port = $PORT\n\
enableCORS = false\n\
headless = true\n\
\n\
" > ~/.streamlit/config.toml