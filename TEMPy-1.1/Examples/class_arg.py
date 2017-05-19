import argparse

class TempyParser(object):

    
    def __init__(self,
                 args=None):
        self.args = args

        
    def generate_args(self):
        parser = argparse.ArgumentParser()
        #
        job_title = parser.add_argument_group()
        job_title.add_argument(
            '-j',
            '--job_title',
            help='Short description of job',
            metavar='Job title',
            type=str,
            dest='job_title',
            default=None,
            required=False)
       #
        inp_map = parser.add_argument_group()
        inp_map.add_argument(
            '-m',
            '--map',
            help=('Input map'),
            metavar='Input map',
            type=str,
            dest='inp_map',
            default=None,
            required=False)
        inp_map.add_argument(
            '-m1',
            '--map1',
            help=('Input map1'),
            metavar='Input map1',
            type=str,
            dest='inp_map1',
            default=None,
            required=False)
	inp_map.add_argument(
            '-m2',
            '--map2',
            help=('Input map2'),
            metavar='Input map2',
            type=str,
            dest='inp_map2',
            default=None,
            required=False)

	search_model = parser.add_argument_group()
        search_model.add_argument(
	    '-p',
            '--pdb',
            help=('Input model'),
            type=str,
            default=None,
	    dest='pdb',
	    required=False)
	map_res = parser.add_argument_group()
        map_res.add_argument(
	    '-r',
            '-res',
            '--map_res',
            help=('Resolution of the map'),
            type=float,
            dest='res',
            default=None)
	map_res.add_argument(
            '-r1',
            '-res1',
            '--map_res1',
            help=('Resolution of the map1'),
            type=float,
            dest='res1',
            default=None)
	map_res.add_argument(
            '-r2',
            '-res2',
            '--map_res2',
            help=('Resolution of the map2'),
            type=float,
            dest='res2',
            default=None)
	map_threshold = parser.add_argument_group()
        map_threshold.add_argument(
            '-t',
            '--thr',
            help='Map threshold value',
            type=float,
            dest='thr',
            default=None)
	map_threshold.add_argument(
            '-t1',
            '--thr1',
            help='Map1 threshold value',
            type=float,
            dest='thr1',
            default=None)
	map_threshold.add_argument(
            '-t2',
            '--thr2',
            help='Map2 threshold value',
            type=float,
            dest='thr2',
            default=None)
	map_threshold.add_argument(
            '-lT',
            '--tlow',
            help='Lower map threshold value',
            type=float,
            default=None)
	map_threshold.add_argument(
            '-hT',
            '--thigh',
            help='Higher map threshold value',
            type=float,
            default=None)
	num_hits = parser.add_argument_group()
        num_hits.add_argument(
            '-hits',
            '--nhits',
            help=('Number of solutions required (1-10000)'),
            type=int,
            dest='num_hits',
            default=100)
	list_inp = parser.add_argument_group()
        list_inp.add_argument(
            '-l',
            '--lfile',
            help=('Input list file'),
            type=str,
            dest='lfile',
            default=100)
	mode = parser.add_argument_group()
        mode.add_argument(
            '-ref',
            '--mode',
            help=('Input mode'),
            type=int,
            dest='mode',
            default=None)
	apix = parser.add_argument_group()
        apix.add_argument(
            '-s',
            '--apix',
            help=('Input grid spacing'),
            type=float,
            dest='apix',
            default=None)
	flag_filt = parser.add_argument_group()
        flag_filt.add_argument(
            '-filt',
            '--filt',
            help=('filter map?'),
            type=bool,
            dest='filt',
            default=True)
	flag_scale = parser.add_argument_group()
        flag_scale.add_argument(
            '-scale',
            '--scale',
            help=('Scale amplitudes?'),
            type=bool,
            dest='scale',
            default=True)
        #
        self.args = parser.parse_args()

