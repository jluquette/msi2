# Anything related to running the Cosmos pipeline on orchestra.

def make_drmaa_lsf_spec(project_name):
    def drmaa_lsf_spec(jobAttempt):
        """Define the flags to be passed to bsub when jobs are queued."""
        task = jobAttempt.task
    
        if task.time_requirement <= 10:
            queue = "mini"
            max_time = '0:10'
        elif task.time_requirement <= 12*60:
            queue = "short"
            max_time = '12:00'
        else:
            queue = "long"
            max_time = '120:00'

        # Override the previous queue.  mcore allows job times up to 1 month
        if task.cpu_requirement > 1:
            queue = "mcore"
    
        short_name = ''.join(c for c in task.stage.name if c.isalnum())
        lsf_name = "/park/msi/%s/%s" % (project_name, short_name)
    
        # Optional list of hosts to exclude from LSF scheduling.
        # XXX: move exclude_hosts' definition to settings.
        exclude_hosts = [ 'clarinet001-253' ]
        exclusion_string = "select[%s]" % ' && '.join('hname != ' + h for h in exclude_hosts) \
            if len(exclude_hosts) > 0 else ''
    
        # task.mem_req MUST be in MB.  When LSF schedules, it multiplies the
        # memory requirement by the -n parameter.  An extra GB of RAM just in
        # case.
        mem_req = int(task.memory_requirement/task.cpu_requirement) + 1024

        param_str = "-q %s -W %s -R 'rusage[mem=%s] %s' -n %s -J %s" % \
            (queue, max_time, mem_req,
            exclusion_string, task.cpu_requirement, lsf_name)
        return param_str

    return drmaa_lsf_spec
