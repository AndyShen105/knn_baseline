#!/usr/bin/python
# -*- coding: UTF-8 -*-
import threading
import os
import time
import logging



logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    )

setting = [

    (247753, 50, 100, 20, 33670, 0, 0, "/MovieLens/q2.txt", "/MovieLens/q2.txt"),
    (247753, 50, 100, 20, 33670, 0, 0, "/MovieLens/q2.txt", "/MovieLens/q2.txt"),
    (247753, 50, 100, 20, 33670, 0, 0, "/MovieLens/q2.txt", "/MovieLens/q2.txt"),

]
def execute(n_rows, n_col, n_neighbour, n_query, n_sample_query, n_sample_data, data_file, query_file):
    """
    :param cmd:
    :return:
    :rtype: str
    """
    logging.info("run config: n_rows:%d, n_col:%d, n_neighbour:%d, n_query:%d, n_sample_query:%s, n_sample_data:%d, data_file:%s, query_file%s" % (
        n_rows,
        n_col,
        n_neighbour,
        n_query,
        n_sample_query,
        n_sample_data,
        data_file,
        query_file
    ))
    cmd = "./baseline %d %d %d %d %d %d %s %s" % (n_rows,
                                                n_col,
                                                n_neighbour,
                                                n_query,
                                                n_sample_query,
                                                n_sample_data,
                                                data_file,
                                                query_file)

    logging.info("run command: %s" % cmd)
    os.system(cmd+" >> running_log.out")


def run():
    for config in setting:
        n_rows = config[0]
        n_col = config[1]
        n_neighbour = config[2]
        n_query = config[1]
        n_sample_query = config[3]
        n_sample_data = config[4]
        data_file = config[5]
        query_file = config[6]
        threads = []
        t1 = threading.Thread(target=execute,args=(n_rows, n_col, n_neighbour, n_query, n_sample_query, n_sample_data, data_file, query_file))
        threads.append(t1)
        t1.setDaemon(True)
        t1.start()
        time.sleep(10)

if __name__=="__main__":
    run()