import sqlite3
import pandas as pd

# Connect to your database
conn = sqlite3.connect("svdb.db")

# See what tables are inside
tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table';", conn)
print(tables)

# Once you know a table name, you can view it:
df = pd.read_sql_query("SELECT * FROM your_table_name LIMIT 10;", conn)
print(df.head())

# Close when you're done
conn.close()
