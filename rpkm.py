#import norm from bioinfokit
import bioinfokit
from bioinfokit.analys import norm, get_data
#import panda for read and write
import pandas as pd
#import your data
df = pd.read_table('/media/scbb/data6/Abhijit/count_final.txt')
#check your data
df.head(2)
# make gene column as index column
df = df.set_index('Gene')
df.head(2)
#now, normalize raw counts using rpkm method
# gene length must be in bp
nm = norm()
#get rpkm normalized dataframe
nm.rpkm(df=df, gl='length')
rpkm_df = nm.rpkm_norm
rpkm_df.head(2)
#export your data
rpkm_df.to_csv(r'/media/scbb/data6/Abhijit/rpkm_data.csv',index=True, header=True,encoding = 'utf-8',  sep ='\t')
