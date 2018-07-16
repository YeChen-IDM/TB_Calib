from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from dtk.utils.builders.sweep import GenericSweepBuilder
from simtools.SetupParser import SetupParser
from tb.TBCustomReports import *
from tb.add_tbhiv_treat import add_tbhiv_treat
from tb.add_ActiveDiagnostic import add_ActiveDiagnostic
from tb.add_HIVIncidence import add_HIVIncidence
from tb.add_DiagnosticTreatNeg import add_DiagnosticTreatNeg
from tb.add_SimpleHealthSeeking import add_SimpleHealthSeeking
from tb.add_Ramp_DiagnosticTreatNeg import add_Ramp_DiagnosticTreatNeg
from tb.add_cd4diagnostic import add_cd4diagnostic
from tb.add_art import add_art
from tb.add_ResistanceDiagnostic import add_ResistanceDiagnostic
from tb.add_tb_drug_type import add_tb_drug_type
from tb.add_tbhiv_outbreak import add_tbhiv_outbreak
from tb.add_simplehivdiagnostic import add_simplehivdiagnostic

from GujCalibSite import GujCalibSite
from IndiaCalibSite import IndiaCalibSite
# from GujCalibSite_noage import GujCalibSite
from calibtool.CalibManager import CalibManager
from calibtool.Prior import MultiVariatePrior
from calibtool.algorithms.IMIS import IMIS
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from calibtool.plotters.SiteDataPlotter import SiteDataPlotter
from random import randint

sites = [GujCalibSite()]
# sites = [IndiaCalibSite()]



prior = MultiVariatePrior.by_range(
    Base_Infectivity=('linear', 0.01, 0.0329),
    Post_Infection_Acquisition_Multiplier=('linear', 0, 1.0),
    #Immunity_Acquisition_Factor=('linear', 0, 1.0),
    #Seek_Time_Slow= ('linear',low_seek_time_min, low_seek_time_max),
    #Seek_Time_Fast= ('linear', high_seek_time_min, high_seek_time_max)
    #TB_Active_Presymptomatic_Infectivity_Multiplier= ('linear',1e-9, 0.3),
    #TB_Presymptomatic_Rate= ('linear', 0.0014, 0.0333 )



    #TB_Fast_Progressor_Rate= ('linear', 1e-4, 0.2)
)

plotters = [
    #LikelihoodPlotter(combine_sites=True)
    #SiteDataPlotter(combine_sites=True)
]


def sample_point_fn(cb, sample_dimension_values):
    """
    A simple example function that takes a list of sample-point values
    and sets parameters accordingly using the sample-dimension names from the prior.
    Note that more complicated logic, e.g. setting campaign event coverage or habitat abundance by species,
    can be encoded in a similar fashion using custom functions rather than the generic "set_param" or "update_params".
    """

    # TODO: reconcile variable names with Pull Request #687/#733, i.e. function accepts one row of sample_point_table?
    sample_point = prior.to_dict(sample_dimension_values)  # aligns names and values; rounds integer-range_type params

    params_to_update = dict()

    #for sample_dimension_name, sample_dimension_value in sample_point.items():
        # Apply specific logic to convert sample-point dimensions into simulation configuration parameters
        #if sample_dimension_name == 'Seek_Time_Slow':
        #    LowSeek(cb,1.0/sample_dimension_value)
        #    params_to_update[sample_dimension_name] = sample_dimension_value
        #elif sample_dimension_name == 'Seek_Time_Fast':
        #    HighSeek(cb,1.0/sample_dimension_value)
        #    params_to_update[sample_dimension_name] = sample_dimension_value
        #else:
        #    if '_LOG' in sample_dimension_name:
        #        param_name = sample_dimension_name.replace('_LOG', '')
        #        params_to_update[param_name] = pow(10, sample_dimension_value)
        #    else:
        #        params_to_update[sample_dimension_name] = sample_dimension_value

    return cb.update_params(params_to_update)


# Set the default configuration block to HPC (we will run on COMPS)
SetupParser.default_block = 'HPC'

# Choose a name for our experiment
exp_name = 'Gujarat cal short'
# Create a default ConfigBuilder
cb = DTKConfigBuilder.from_defaults('TBHIV_SIM')
#burn in parameters
burn_initial = 1*365
burn_predots = 300*365
#burn_predots = 30*365


dots_start = burn_initial + burn_predots
genexpert_introduction = dots_start + 9*365
length_gene_xpert_ramp = 5*365
To_end_from_DOTS = 35*365
#To_end_from_DOTS = 5*365

#passive seeking rates
high_seek = 1.0/(3*30)
low_seek = 1.0/(7*30)

sens_smear_neg_pre_GL = 0.4
sens_smear_pos_pre_GL  = 0.7

sens_smear_neg_pre_GH = 0.4
sens_smear_pos_pre_GH  = 0.7


sens_smear_neg_GXL = 0.7
sens_smear_pos_GXL = 0.7
sens_smear_pos_GXH  = 0.99
sens_smear_neg_GXH = 0.9

sens_resistance_L = 0.5
sens_resistance_H = 0.7
specificity_resistance = 0.95

# set potentially calibratable param
#cb.set_param('TB_Active_Presymptomatic_Infectivity_Multiplier', 1e-9)
#cb.set_param('TB_Presymptomatic_Rate', 1/(30*3))

cb.set_param('TB_Active_Presymptomatic_Infectivity_Multiplier', 0.3)
cb.set_param('TB_Presymptomatic_Rate', 1/(30*6))

# set underlying parameters
cb.set_param('Listed_Events', ['Blackout', 'TBTestPreDOTSLow','TBTestPreDOTSHigh','TBDS_Positive',
                               'TBTestDOTSHigh', 'TBTestDOTSLow', 'Active_Smear_Neg'])
cb.set_param('Simulation_Duration', burn_initial + burn_predots + To_end_from_DOTS)
cb.set_param('Base_Population_Scale_Factor', 20000)
#cb.set_param('Base_Infectivity', 0.0219)  # set as calibration parameter
cb.set_param('Demographics_Filenames', ['Base_Demog_India.json', 'Base_Overlay_India.json'])
#cb.set_param('x_Other_Mortality', 0)
cb.set_param('x_Birth', 1.0)

# dynamics during initial burn in (natural dynamics only)
#initial seeding loaded in from_defaults


# add an initial outbreak
add_tbhiv_outbreak(cb, 0.02, 'TB')

# second burn in phase PREDOTS treatment and diagnostics


add_DiagnosticTreatNeg(cb, ['TBTestDOTSHigh'], sens_smear_pos_pre_GH, sens_smear_neg_pre_GH, treatment_fraction=0.8,
                       start_day=dots_start, duration=genexpert_introduction - dots_start, property_restrictions_list=['Care_Quality:High'])
add_DiagnosticTreatNeg(cb, ['TBTestDOTSLow'], sens_smear_pos_pre_GL, sens_smear_neg_pre_GL, treatment_fraction=0.8,
                       start_day=dots_start, duration=genexpert_introduction-dots_start, property_restrictions_list=['Care_Quality:Low'])

add_Ramp_DiagnosticTreatNeg(cb,['TBTestDOTSHigh'],length_gene_xpert_ramp, sens_smear_pos_GXH, sens_smear_neg_GXH,
                            sens_smear_pos_pre_GH, sens_smear_neg_pre_GH, 0.8, treatment_fraction=0.8,
                            pos_event= 'ProviderOrdersTBTest',
                            start_day= genexpert_introduction, duration= -1, property_restrictions_list= ['Care_Quality:High'] )

add_Ramp_DiagnosticTreatNeg(cb, ['TBTestDOTSLow'], length_gene_xpert_ramp, sens_smear_pos_GXL, sens_smear_neg_GXL,
                            sens_smear_pos_pre_GL, sens_smear_neg_pre_GL, 0.8, treatment_fraction=0.8,
                            pos_event= 'ProviderOrdersTBTest',
                            start_day=genexpert_introduction, duration=-1,
                            property_restrictions_list=['Care_Quality:Low'])

add_ResistanceDiagnostic(cb, ['ProviderOrdersTBTest'],sens_resistance_L, specificity_resistance, neg_event = 'TBDS_Positive',
                         treatment_fraction= 0.5, treatment_fraction_negative_test = 1.0,
                         start_day= genexpert_introduction, property_restrictions_list=['Care_Quality:Low'] )

add_ResistanceDiagnostic(cb, ['ProviderOrdersTBTest'],sens_resistance_H, specificity_resistance, neg_event = 'TBDS_Positive',
                         treatment_fraction= 0.5, treatment_fraction_negative_test = 1.0,
                         start_day= genexpert_introduction, property_restrictions_list=['Care_Quality:High'] )

add_tbhiv_treat(cb, 'DOTSHQ', ['TBTestPositive','TBDS_Positive' ], start_day= dots_start,
                latent_multiplier=0, property_restrictions_list=['Care_Quality:Low'])

add_tbhiv_treat(cb, 'DOTSLQ', ['TBTestPositive', 'TBDS_Positive'], start_day= dots_start,
                latent_multiplier=0, property_restrictions_list=['Care_Quality:High'])

add_tbhiv_treat(cb, 'DOTSMDR', ['TBMDRTestPositive'], start_day=dots_start,
                latent_multiplier=0 )
#ART introduction in 2007
#maybe include ART dropout at a rate

#add in appropriate DLLs
# add_tb_report(cb, additional_events= ['Active_Smear_Neg'],stop_year= 2000, min_age_yrs=0, max_age_yrs= 200, type= 'Report_TBHIV_ByAge' )
add_tb_report(cb, stop_year= 2000, min_age_yrs=0, max_age_yrs= 200, type= 'Report_TBHIV_ByAge' )
#add_tb_drug_type(cb, 'PreDOTSHigh',180,0.5,0.03,0,0.10,0.02)
#add_tb_drug_type(cb, 'PreDOTSLow',180,0.5,0.03,0,0.10,0.02)
add_tb_drug_type(cb, 'DOTSHQ',180,0.8,0.03,0,0.10,0.02)
add_tb_drug_type(cb, 'DOTSLQ',180,0.5,0.03,0,0.10,0.02)
add_tb_drug_type(cb, 'DOTSMDR',180,0.5,0.03,0,0.10,0.02)

test = cb.get_param('Listed_Events')
cb.set_param('logLevel_default', 'ERROR')
cb.disable('Default_Reporting')
#cb.set_param('Run_Number', randint(0, 1000) )


#next_point_kwargs = dict(initial_samples=5000,
#                         samples_per_iteration=500,
#                         n_resamples=1000)

cb.set_param('Run_Number', randint(0, 10) )


next_point_kwargs = dict(initial_samples=4,
                         samples_per_iteration=3,
                         n_resamples=4)


calib_manager = CalibManager(name='Gujarat_short',
                             config_builder=cb,
                             map_sample_to_model_input_fn=sample_point_fn,
                             sites=sites,
                             next_point=IMIS(prior, **next_point_kwargs),
                             sim_runs_per_param_set=1,
                             max_iterations=2,
                             plotters=plotters)
run_calib_args = {}

if __name__ == "__main__":
    SetupParser.init(selected_block=SetupParser.default_block)
    calib_manager.run_calibration()
