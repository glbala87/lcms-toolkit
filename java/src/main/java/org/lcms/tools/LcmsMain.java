package org.lcms.tools;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Main entry point for LCMS command-line tools.
 */
@Command(
    name = "lcms",
    mixinStandardHelpOptions = true,
    version = "LCMS Tools 1.0.0",
    description = "LC-MS data analysis command-line tools",
    subcommands = {
        InfoCommand.class,
        ConvertCommand.class,
        PeakPickCommand.class,
        XicCommand.class,
        CommandLine.HelpCommand.class
    }
)
public class LcmsMain implements Runnable {

    @Option(names = {"-v", "--verbose"}, description = "Enable verbose output")
    private boolean verbose;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new LcmsMain()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public void run() {
        // When called without subcommand, show help
        CommandLine.usage(this, System.out);
    }

    public boolean isVerbose() {
        return verbose;
    }
}
