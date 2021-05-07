import os




def simplify_path(path,remove_gz=True):
    """Removes dir and extension from a filepath.
        checks if file has an e
    """
    name,ext= os.path.splitext(os.path.basename(path))

    if remove_gz & (ext=='.gz'):
        name=  os.path.splitext(name)[0]

    return name

def cat_files(files,outfilename,gzip=False):
    """ cat files in python
        gzip: compress outfile
        set to false when cat files that are already gzipped.
    """

    import shutil

    if gzip:
        import gzip as gz
        outhandle= gz.open
    else:
        outhandle = open

    with outhandle(outfilename, 'wb') as f_out:
        for f in files:
            with open(f, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

def convert_percentages(df):
    """Convet all columns with strings and % at the end to percentages
    """
    for col in df.columns:
        if df.dtypes[col]=='object':
            if df[col].iloc[0].endswith('%'):
                df.loc[:,col]= df[col].str.rstrip('%').astype('float') / 100.0


def symlink_relative(files,input_dir,output_dir):
    """create symlink with and adjust for relative path"""


    input_dir_rel= os.path.relpath(input_dir, output_dir)

    for f in files:
        os.symlink(os.path.join(input_dir_rel,f),
                   os.path.join(output_dir,f))

def pandas_concat(input_tables,output_table,sep='\t',index_col=0,axis=0,
                  read_arguments=None,save_arguments=None,concat_arguments=None):
    """
        Uses pandas to read,concatenate and save tables using pandas.concat
    """

    import pandas as pd

    if read_arguments is None:
        read_arguments={}
    if save_arguments is None:
        save_arguments={}
    if concat_arguments is None:
        concat_arguments={}

    if type(input_tables) == str:
        input_tables= [input_tables]

    Tables= [pd.read_csv(file,index_col=index_col,sep=sep,**read_arguments) for file in input_tables]

    out= pd.concat(Tables,axis=axis,**concat_arguments).sort_index()

    del Tables

    out.to_csv(output_table,sep=sep,**save_arguments)

    del out
