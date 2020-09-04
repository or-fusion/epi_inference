import sys
import pandas as pd
import csv

df = pd.read_csv(sys.argv[1], encoding="ISO-8859-1", dtype={'FIPS':str})
print(df)

df = df.rename(columns={'dates':'Date'})
print(df)

df = pd.pivot_table(df, index='Date', columns='FIPS', values='transmissions')
print(df)

df.to_csv('foo.csv', quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")
