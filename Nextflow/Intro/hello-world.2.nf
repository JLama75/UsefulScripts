#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to standard out
 */
params.greeting= "Bonjour le monde!"

process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output: 
        //stdout
        path "${greeting}-output.txt"
    
    """
    echo '$greeting' > '$greeting-output.txt'
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
    greeting_ch = Channel.of('Hello','Bonjour','Holà')
    // create a channel for inputs from a file
    //greeting_ch = Channel.fromPath(params.input_file).splitText() { it.trim() }
    //{ it.trim() } trim the white line for each split text
    // emit a greeting. within the curly braces is the groovy code/syntax
    // emit a greeting
    sayHello(greeting_ch)

    //convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
// nextflow run hello-world.nf -resume --greeting 'Hola!'
// nextflow run hello-world.nf -ansi-log false
