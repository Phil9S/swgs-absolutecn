import argparse
import re
import ruamel.yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString, DoubleQuotedScalarString

yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.explicit_start = True

def get_input(x,y):
    if y == 'bins':
        bin_input = input(f'Set {y} (current: {x}): ')
        if bin_input:
            new_input = list(map(int,bin_input.strip().split(',')))
        else:
            return x
    else:
        new_input = (input(f'Set {y} (current: {x}): '))
    if not new_input:
        return x
    else:
        return new_input

parser = argparse.ArgumentParser(description='Update pipeline config/job submission yamls')
parser.add_argument("-c","--config", \
    help="config to update", \
    type=str, \
    required=True, \
    nargs=1, \
    choices=['config','slurm','pbs','lsf','local'])
parser.add_argument("--advanced", help="Modify advanced cluster configuration",
                    action="store_true")
args = parser.parse_args()

if args.config[0] == "config":
    # Config file updates
    print('Update config.yaml values (Return no value to keep the current value)')
    file_name = 'config/config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    # Update config values
    config['samplesheet'] = DoubleQuotedScalarString(get_input(config['samplesheet'],'samplesheet'))
    config['bins'] = list(get_input(config['bins'],'bins'))
    config['out_dir'] = DoubleQuotedScalarString(get_input(config['out_dir'],'out_dir'))
    config['project_name'] = DoubleQuotedScalarString(get_input(config['project_name'],'project_name'))
    config['af_cutoff'] = float(get_input(config['af_cutoff'],'af_cutoff'))
    config['use_seed'] = DoubleQuotedScalarString(get_input(config['use_seed'],'use_seed'))
    config['seed_val'] = DoubleQuotedScalarString(get_input(config['seed_val'],'seed_val'))
    config['filter_underpowered'] = DoubleQuotedScalarString(get_input(config['filter_underpowered'],'filter_underpowered'))
    # Save updated yaml
    with open('config/config.yaml', 'w') as fp:
        yaml.dump(config, fp)
elif args.config[0] == "slurm":
    # Slurm profile updates
    print('Update slurm profile values (Return no value to keep the current value)')
    file_name = 'profile/slurm/config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    # Update jobs
    config['jobs'] = int(get_input(config['jobs'],'jobs'))
    # Save updated yaml
    with open('profile/slurm/config.yaml', 'w') as fp:
        yaml.dump(config, fp)
    # Cluster config
    file_name = 'profile/slurm/cluster_config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    # Update slurm cluster config
    config['__default__']['account'] = get_input(config['__default__']['account'],'account')
    config['__default__']['partition'] = get_input(config['__default__']['partition'],'partition')
    config['__default__']['time'] = int(get_input(config['__default__']['time'],'default time'))
    config['__default__']['mem'] = int(get_input(config['__default__']['mem'],'default mem'))
    # Advanced options
    if args.advanced:
        config['__default__']['nodes'] = int(get_input(config['__default__']['nodes'],'default nodes'))
        config['__default__']['ntasks'] = int(get_input(config['__default__']['ntasks'],'default ntasks'))
        config['__default__']['cpus-per-task'] = int(get_input(config['__default__']['cpus-per-task'],'default threads'))
        config['__default__']['output'] = get_input(config['__default__']['output'],'output log')
        config['__default__']['error'] = get_input(config['__default__']['error'],'error log')
    # High resource jobs
    config['rel_to_abs']['time'] = int(get_input(config['rel_to_abs']['time'],'high-resource time'))
    config['gridsearch_filter']['time'] = config['rel_to_abs']['time']
    config['rel_to_abs']['mem'] = int(get_input(config['rel_to_abs']['mem'],'high-resource mem'))
    config['gridsearch_filter']['mem'] = config['rel_to_abs']['mem']
    config['rel_to_abs']['cpus-per-task'] = int(get_input(config['rel_to_abs']['cpus-per-task'],'high-resource threads'))
    config['gridsearch_filter']['cpus-per-task'] = config['rel_to_abs']['cpus-per-task']
    if args.advanced:
        config['rel_to_abs']['nodes'] = int(get_input(config['rel_to_abs']['nodes'],'high-resource nodes'))
        config['gridsearch_filter']['nodes'] = config['rel_to_abs']['nodes']
        config['rel_to_abs']['ntasks'] = int(get_input(config['rel_to_abs']['ntasks'],'high-resource ntasks'))
        config['gridsearch_filter']['ntasks'] = config['rel_to_abs']['ntasks']
    # Save updated yaml
    with open('profile/slurm/cluster_config.yaml', 'w') as fp:
        yaml.dump(config, fp)
elif args.config[0] == "local":
    print('Update local profile values (Return no value to keep the current value)')
    file_name = 'profile/local/config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    config['jobs'] = int(get_input(config['jobs'],'jobs'))
    # Save updated yaml
    with open('profile/local/config.yaml', 'w') as fp:
        yaml.dump(config, fp)
    file_name = 'profile/local/cluster_config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    config['cpus-per-task'] = int(get_input(config['cpus-per-task'],'default threads'))
    # Save updated yaml
    with open('profile/local/cluster_config.yaml', 'w') as fp:
        yaml.dump(config, fp)
elif args.config[0] == "pbs":
    print('Update pbs-torque profile values (Return no value to keep the current value)')
    file_name = 'profile/pbs/config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    config['jobs'] = int(get_input(config['jobs'],'jobs'))
    # Save updated yaml
    with open('profile/pbs/config.yaml', 'w') as fp:
        yaml.dump(config, fp)
    # default resources
    file_name = 'profile/pbs/cluster_config.yaml'
    config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
    config['__default__']['A'] = get_input(config['__default__']['A'],'account/user')
    #config['__default__']['q'] = get_input(config['__default__']['q'],'partition')
    # default resources
    resourceSplit = config['__default__']['l'].split(',')
    if args.advanced:
        proc = resourceSplit[0].split(':')
        proc_n = proc[0].split('=')
        proc_n[1] = get_input(proc_n[1],'default nodes')
        proc[0] = '='.join([str(elem) for elem in proc_n])
        proc_t = proc[1].split('=')
        proc_t[1] = get_input(proc_t[1],'default threads')
        proc[1] = '='.join([str(elem) for elem in proc_t])    
        resourceSplit[0] = ':'.join([str(elem) for elem in proc])

    time = resourceSplit[1].split('=')
    time[1] = get_input(time[1],'default wallclock')
    resourceSplit[1] = '='.join([str(elem) for elem in time])
    
    mb = resourceSplit[2].split('=')
    mb[1] = get_input(mb[1],'default mem')
    resourceSplit[2] = '='.join([str(elem) for elem in mb])
    config['__default__']['l'] = ','.join([str(elem) for elem in resourceSplit])
    # high resources
    resourceSplit = config['rel_to_abs']['l'].split(',')
    if args.advanced:    
        proc = resourceSplit[0].split(':')
        proc_n = proc[0].split('=')
        proc_n[1] = get_input(proc_n[1],'high-resources nodes')
        proc[0] = '='.join([str(elem) for elem in proc_n])

        proc_t = proc[1].split('=')
        proc_t[1] = get_input(proc_t[1],'high-resource threads')
        proc[1] = '='.join([str(elem) for elem in proc_t])
        resourceSplit[0] = ':'.join([str(elem) for elem in proc])
    
    time = resourceSplit[1].split('=')
    time[1] = get_input(time[1],'high-resource wallclock')
    resourceSplit[1] = '='.join([str(elem) for elem in time])

    mb = resourceSplit[2].split('=')
    mb[1] = get_input(mb[1],'high-resource mem')
    resourceSplit[2] = '='.join([str(elem) for elem in mb])
    config['rel_to_abs']['l'] = ','.join([str(elem) for elem in resourceSplit])
    config['gridsearch_filter']['l'] = config['rel_to_abs']['l']
    # Save updated yaml
    with open('profile/pbs/cluster_config.yaml', 'w') as fp:
        yaml.dump(config, fp)
elif args.config[0] == "lsf":
    print("lsf is not currenty implemented!")

