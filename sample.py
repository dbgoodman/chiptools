import string
import re
from collections import namedtuple

def num_to_well(num):
    num = int(num)
    col = (num-1) / 8 + 1
    row = (num-1) % 8
    row = string.ascii_uppercase[row]
    return '%s%0.2d' % (row,col)
    
def well_to_num(well):
    well = well.upper()
    row,col = re.match(r'(\w)(\d+)',well).groups()
    col = int(col)
    return 8*(col-1)+(string.ascii_uppercase.find(row)+1)
   
def parse_sanger_id(sanger_id_string, regex_str=r'([\w-]+)-(\d+)', 
                    **kwargs):
    sanger_pos = namedtuple('sanger_pos', ['plate', 'num', 'well', 'id'])
    (plate,num) = re.match(regex_str,sanger_id_string).groups()
    well = num_to_well(num)
    
    return sanger_pos(plate, num, well, sanger_id_string)

