"""

"""
import multiprocessing as mp
from Queue import Empty
import time

class Worker(mp.Process):
    def __init__(self, results_queue, finished_queue, ID, **opts):
        super(Worker, self).__init__(**opts)
        self.ID = ID
        self.q = results_queue
        self.fq = finished_queue
    def run(self):
        print "RUNNING %s" % self.ID
        res = self._target()
        self.q.put((self.ID, res))
        print "FINISHED %s" % self.ID
        self.fq.put(self.ID)
        return

class Jobs:
    def __init__(self, nr_procs=2, wait_time=0.05):
        self.nr_procs = nr_procs
        self.jobs = []
        self.results_queue = mp.Queue()
        self.finished_queue = mp.Queue()
        self.nr_jobs = 0
        self.wait_time = wait_time
    def job(self, f):
        self.jobs.append( Worker(self.results_queue, self.finished_queue, self.nr_jobs, target=f, name=f.func_name) )
        self.nr_jobs += 1
    def count_running(self):
        nr_running = 0
        for i,p in enumerate(self.jobs):
#            print "Is worker %s alive? %s" % (i,p.is_alive())
            if p.is_alive():
                nr_running += 1
        nr_running -= self.results_queue.qsize()
        return nr_running
    def run_parallel(self):
        for p in self.jobs:

            while True:
                capacity = self.nr_procs - self.count_running()
#                print "Capacity = %s" % capacity
                if capacity > 0:
                    p.start()
                    break
                else:
                    finished_IDs = []
                    while True:
                        try:
                            ID = self.finished_queue.get(block=False)
                            finished_IDs.append(ID)
                        except Empty:
                            break
                    """
                    print "finished IDs = %s" % finished_IDs
                    for pf in self.jobs:
                        if pf.ID in finished_IDs:
                            print "JOIN %s" % pf.ID
                            pf.join(self.wait_time)
                            print "JOINED %s" % pf.ID
                    """
                    time.sleep(self.wait_time)
#            p.start()
        for p in self.jobs:
            p.join()
            assert p.exitcode == 0, "some process failed"
        results = []
        while True:
            results.append(self.queue.get(False))
            if self.queue.empty():
                break
        # sort results in the same order as the jobs were submitted
        results.sort(key=lambda tup: tup[0])
        results = [r[1] for r in results]
        # remove all jobs
        self.jobs = []
        return results
