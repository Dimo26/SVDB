import sqlite3
import pandas as pd

conn = sqlite3.connect("test_db.db")
cursor = conn.cursor()

# Check what tables exist
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
print("Tables:", tables)

# Query the SVDB table (not test_db)
df = pd.read_sql_query("SELECT * FROM SVDB LIMIT 10;", conn)
print("\nFirst 10 rows from SVDB table:")
print(df.head(10))

# Check total number of variants
cursor.execute("SELECT COUNT(*) FROM SVDB;")
count = cursor.fetchone()[0]
print(f"\nTotal variants in database: {count}")

# Show unique variant types
cursor.execute("SELECT DISTINCT var FROM SVDB;")
var_types = cursor.fetchall()
print(f"\nVariant types: {[v[0] for v in var_types]}")



# Count variants by type

conn.close()