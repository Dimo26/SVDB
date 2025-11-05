import sqlite3
import pandas as pd

conn = sqlite3.connect("SVDB.db")
cursor = conn.cursor()

cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
print("Tables:", tables)

df = pd.read_sql_query("SELECT * FROM SVDB LIMIT 10;", conn)
print(df.head())

conn.close()
