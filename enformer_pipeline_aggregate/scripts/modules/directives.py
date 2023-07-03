

def parsl_directives(use_parsl, parsl_parameters):

    import parsl

    if use_parsl == True:
        if parsl_parameters['provider'] == 'highthroughput':
            if parsl_parameters['hpc'] == 'polaris':
                print(f"INFO - Using {parsl_parameters['provider']} parsl configuration for polaris: {use_parsl}")
                import parslConfiguration
                parsl.load(parslConfiguration.polaris_htParslConfig(params=parsl_parameters))
            elif parsl_parameters['hpc'] == 'theta':
                print(f'INFO - Using parsl {parsl_parameters["provider"]} configuration for theta: {use_parsl}')
                import parslConfiguration
                parsl.load(parslConfiguration.theta_htParslConfig(params=parsl_parameters))
            elif parsl_parameters['hpc'] == 'beagle3':
                print(f'INFO - Using parsl {parsl_parameters["provider"]} configuration for beagle3: {use_parsl}')
                import parslConfiguration
                parsl.load(parslConfiguration.beagle3_htParslConfig(params=parsl_parameters))
        elif parsl_parameters['provider'] == 'local':
            if parsl_parameters['hpc'] == 'polaris':
                print(f'INFO - Using parsl {parsl_parameters["provider"]} configuration for polaris: {use_parsl}')
                import parslConfiguration
                parsl.load(parslConfiguration.polaris_localParslConfig(params=parsl_parameters))
            elif parsl_parameters['hpc'] == 'theta':
                print(f'INFO - Using parsl {parsl_parameters["provider"]} configuration for theta: {use_parsl}')
                import parslConfiguration
                parsl.load(parslConfiguration.theta_localParslConfig(params=parsl_parameters))
            elif parsl_parameters['hpc'] == 'beagle3':
                print(f'INFO - Using parsl {parsl_parameters["provider"]} configuration for beagle3: {use_parsl}')
                import parslConfiguration
                parsl.load(parslConfiguration.beagle3_localParslConfig(params=parsl_parameters))

    elif use_parsl == False:
        return(0)
        # if parsl_parameters['hpc'] == 'polaris':
        #     print(f'INFO - Using parsl configuration for polaris: {use_parsl}')
        #     import parslConfiguration
        #     parsl.load(parslConfiguration.polaris_localParslConfig(params=parsl_parameters))
        # elif parsl_parameters['hpc'] == 'theta':
        #     print(f'INFO - Using parsl configuration for theta: {use_parsl}')
        #     import parslConfiguration
        #     parsl.load(parslConfiguration.theta_localParslConfig(params=parsl_parameters))
        
    return(0)