import sqlite3
import json

def export_view_data_to_json(db_path, view_name, json_file_path):
    try:
        # Connect to the SQLite database
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Query the view to get its contents
        cursor.execute(f"SELECT * FROM {view_name}")
        rows = cursor.fetchall()
        
        # Get column names
        column_names = [description[0] for description in cursor.description]
        
        # Convert rows into a list of dictionaries
        data = [dict(zip(column_names, row)) for row in rows]
        
        # Write to a JSON file
        with open(json_file_path, 'w') as json_file:
            json.dump(data, json_file, indent=4)
        
        print(f"Data from view '{view_name}' has been exported to {json_file_path}")
    
    except Exception as e:
        print(f"An error occurred: {e}")
    
    finally:
        if conn:
            conn.close()




# Usage
virus = "Nipah"
db_file = f"OutputData/{virus}/{virus}.db"  # Replace with your database file path
view_name = "tblGBRefs"       # Replace with your view name
output_json = f"OutputData/{virus}/{view_name}.json"   # Desired output JSON file name

export_view_data_to_json(db_file, view_name, output_json)
