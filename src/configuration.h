
    /// Constructor with reasonable default parameters
    general(void) :
      paramFilename("oplsaa.prm"), outputFilename("%i.out"),
      input(input_types::TINKER), output(output_types::TINKER),
      task(config::tasks::SP), energy_interface(interface_types::OPLSAA),
      preopt_interface(interface_types::ILLEGAL),
      verbosity(1U)