import sqlite3 as sql

class KEGGDB(object):
    """docstring for KEGGDB"""
    def __init__(self, db):
        super(KEGGDB, self).__init__()
        self.db = db
        self.connect(db)

    def connect(self, db):
        self.con = sql.connect(db)
        print(db)
        print('DB Connected!')

    def fetchall():
        with con:
            cur = con.cursor()    
            cur.execute("SELECT * FROM Cars")

            rows = cur.fetchall()

            for row in rows:
                print row


if __name__ == '__main__':
  kegg_db = KEGGDB('/Volumes/userdata/Joonyong/bsf-server-gse/KEGGDB_2016-08-30.db')
