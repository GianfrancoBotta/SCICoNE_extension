import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import polars as pl
import re

from sklearn.base import RegressorMixin
from sklearn.linear_model import LinearRegression
from typing import Union, Callable, Optional
   
   
def correct_depth(
    depth: Union[pd.DataFrame, pl.DataFrame],
    corr_vars: list,
    outpath: Optional[str] = None,
) -> Union[pd.DataFrame, pl.DataFrame]:
    """
    Perform different types of correction on depth of genomic data.

    Parameters
    ----------
    depth : pandas.DataFrame | polars.DataFrame
        Dataset containing depth per bin and columns to correct the depth.
    corr_vars : str
        List of variables to correct depth, options are: 'gc_bias', 'bin_length'
    df2 : pandas.DataFrame | polars.DataFrame
        Second dataset, which contains the predictor variable.
    outpath : str
        Directory to save a tsv file of the new dataframe, default to None.
        
    Returns
    -------
    Union[pd.DataFrame, pl.DataFrame]
        depth dataframe corrected with the given variables.
    """
    # Convert pandas to polars
    if not isinstance(depth, pl.DataFrame):
        depth = pl.from_pandas(depth)
    
    # Correct depth
    expr = pl.lit(1)

    if 'gc_bias' in corr_vars and 'gc_bias' in depth.columns:
        expr = expr * pl.col('gc_bias')
    if 'bin_length' in corr_vars and 'bin_length' in depth.columns:
        expr = expr * pl.col('bin_length')

    depth = depth.with_columns(expr.alias("corr_factor"))
    
    # Add normalized depth to df
    depth = depth.with_columns(
        pl.when(pl.col("corr_factor").is_not_null() & (pl.col("corr_factor") != 0))
        .then(pl.col("depth") / pl.col("corr_factor"))
        .otherwise(0)
        .alias('_'.join(corr_vars + ['norm_depth']))
    )
    
    if outpath is not None:
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        depth.write_csv(outpath, separator='\t')
        
    return(depth)

def gtf2df(gtf: str, 
           engine: str = 'polars',
           columns: list = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
           fields: list = ['gene_id', 'transcript_id', 'gene_type']) -> Union[pd.DataFrame, pl.DataFrame]:
    """
    Convert a GTF file to a pandas/polars DataFrame.
    
    Parameters
    ----------
    gtf : str
        Path to the GTF file.
    engine : str
        String to specify if the user wants a Pandas or a Polars dataframe, options are 'pandas' and 'polars', default to 'polars'.
    columns : list
        Column names of the gtf files, default to ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'].
    fields : list
        Fields to split the last column of the gtf file, default to ['gene_id', 'transcript_id', 'gene_type'].
    Returns
    -------
    Union[pd.DataFrame, pl.DataFrame]
    """

    df = pd.read_csv(gtf, sep='\t', header=None, comment='#')
    df.columns = columns
    
    for field in fields:
        df[field] = df['attribute'].apply(lambda x: re.findall(rf'{field} "([^"]*)"', x)[0] if rf'{field} "' in x else '')

    df = df.infer_objects(copy=False)

    df.drop('attribute', axis=1, inplace=True)

    if(engine == 'polars'):
        df = pl.from_pandas(df)
    
    return df
   
def regress(
    df1: Union[pd.DataFrame, pl.DataFrame],
    var1: str, 
    df2: Union[pd.DataFrame, pl.DataFrame],
    var2: str,
    xtransform: Optional[Callable] = None,
    ytransform: Optional[Callable] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    legend_loc: str = 'best',
    R_loc: list = [0.05, 0.95],
    model: Optional[RegressorMixin] = None,
    plotdir: Optional[str] = None,
) -> RegressorMixin:
    """
    Perform regression between variables from two DataFrames (pandas or polars).

    Parameters
    ----------
    df1 : pandas.DataFrame | polars.DataFrame
        First dataset, which contains the predictor variable.
    var1 : str
        Column name from df1, the predictor variable.
    df2 : pandas.DataFrame | polars.DataFrame
        Second dataset, which contains the predictor variable.
    var2 : str
        Column name from df2, the response variable.
    xtransform : callable, optional
        Function to transform x (e.g. np.log), default to None.
    ytransform : callable, optional
        Function to transform y (e.g. np.log), default to None.
    xlabel : str
        Xlabel for the plot, default to var1.
    ylabel : str
        Ylabel for the plot, default to var2.
    legend_loc : str
        Location of the legend in the plot, default to 'best'.
    R_loc : list
        Location of the R^2 value in the plot, default to upper left [0.05, 0.95].
    model : sklearn.base.RegressorMixin
        Any sklearn-compatible regression model (e.g. LinearRegression, Ridge, Lasso, etc.), default to sklearn.linear_model.LinearRegression().
    plotdir : str
        Directory to save the output, default to None if plots are not needed.
        
    Returns
    -------
    np.sklearn.base.RegressorMixin
        The fitted model on the passed data.
    """
    # Convert pandas to polars
    if not isinstance(df1, pl.DataFrame):
        df1 = pl.from_pandas(df1)
    if not isinstance(df2, pl.DataFrame):
        df2 = pl.from_pandas(df2)
        
    # Set default model
    if model is None:
        model = LinearRegression()
    
    # Fit the linear model
    X = df1.select(var1)
    if xtransform is not None:
        X = xtransform(X)
    y = df2.select(var2)
    if ytransform is not None:
        y = ytransform(y)
    
    model.fit(X, y)
    y_pred = model.predict(X)
    R = model.score(X,y)
    
    # Plot scatter and regression line
    plt.figure()
    plt.scatter(X, y, color='#3182bd', label='Data')
    plt.plot(X, y_pred, color='#de2d26', label='Regression line')
    plt.text(R_loc[0], R_loc[1], f"$R^2 = {R:.3f}$", transform=plt.gca().transAxes,
            fontsize=12, verticalalignment='top')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=legend_loc)
    if plotdir is not None:
        os.makedirs(plotdir, exist_ok=True)
        plt.savefig(
            fname=os.path.join(plotdir, 'gc_pct_vs_norm_depth.png'),
            transparent=False,
            dpi=300,
            format='png'
        )
    plt.close()
    
    return(model) 