import csv
import datetime

def is_number(s):
    """
    Checks if the string 's' can be converted to a float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def ac_logger_reader(file_in, col_to_read=[3, ], power_factor=0.8):
    """
    reads data exported to a csv file from a HIOKI ac data logger
    inputs:
        file_in         path to a csv file to read
        col_to_read     List of columns to read
    return:
        data_dict       dict containing a list of each data type specified in col_to_read
    """

    input_data = csv.reader(open(file_in))
    # skip the header row
    row = input_data.__next__()
    # read in data
    data_dict = {}
    key_list = ['start_time','end_time',]
    for row in input_data:
        if row[1] == 'CH comment':
            for key in row[2:]:
                key_list.append(key + ' ')
            continue
        if row[1] == 'Property':
            for i in range(2, len(key_list)):
                key_list[i] += row[i]
            for key in key_list:
                data_dict[key] = []
            continue
        if row[0] is not '':
            for i in range(0, len(row)):
                if is_number(row[i]):
                    row[i] = float(row[i])
                data_dict[key_list[i]].append(row[i])
    for i in range(0,len(data_dict['start_time'])):
        temp = data_dict['start_time'][i].split(sep='/')
        temp1 = temp[-1].split(sep=' ')
        temp[-1] = temp1[0]
        temp1 = temp1[1].split(':')
        temp.extend(temp1)
        data_dict['start_time'][i] = datetime.datetime(int(temp[2]), int(temp[0]), int(temp[1]), int(temp[3]), int(temp[4]))
    for i in range(0,len(data_dict['end_time'])):
        temp = data_dict['end_time'][i].split(sep=':')
        temp1 = temp[-1].split(' ')
        temp[-1] = temp1[0]
        temp[0] = int(temp[0])
        if temp1[-1] == 'PM' and temp[0]<12:
            temp[0] += 12
        data_dict['end_time'][i] = datetime.time(temp[0],int(temp[1]),int(temp[2]))

    return(data_dict)