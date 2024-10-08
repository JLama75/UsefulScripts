#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to standard out
 */
params.gretting= "Bonjour le monde!"

process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output: 
        //stdout
        path 'output.txt'
    
    """
    echo '$greeting' > output.txt
    """
}

/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}

workflow {

    //create a channel for inputs. Build me a channel that has the element hello world!
    greeting_ch = Channel.of(params.greeting)
    // emit a greeting
    sayHello(greeting_ch)

    //convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
