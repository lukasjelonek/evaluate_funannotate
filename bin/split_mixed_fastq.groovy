#!/usr/bin/env groovy

@Grapes(
  @Grab(group='io.reactivex.rxjava3', module='rxjava', version='3.0.4')
)
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.nio.file.Files
import java.nio.file.Paths
import io.reactivex.rxjava3.core.*
import groovy.transform.CompileStatic

def cli = new CliBuilder(
        usage: "split_mixed_fastq.groovy -f <fastq.gz file> -p <output_prefix>",
        header: '\n Splits a fastq file into multiple files based on the pattern `/\\d+$` in the headers. ' +
                'For each distinct number a new file with the scheme ' +
                '$prefix_$number.fastq.gz is created.\n',
        footer: ''
        )
cli.f(longOpt:'fastq', args:1, argName: 'fastq', 'The fastq.gz file')
cli.p(longOpt:'prefix', args:1, argName: 'prefix', 'The output prefix')
cli.h(longOpt:'help', 'Show help')

def options = cli.parse(args)

if (!options || options.h || !options.f || !options.p) {
  println(cli.usage())
  System.exit(1)
}

String file = options.f
String prefix = options.p
process(file, prefix)

@CompileStatic
void process(String file, String prefix) {
    BufferedReader input = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(Paths.get(file)))))

    Map<String, BufferedWriter> writeFiles = [:]

    String cur_id = null
    input.eachLine{line, number ->
        if ((number-1) % 4 == 0) {
            String[] s = line.split("/")
            cur_id = s[-1]
            line = line.replace(" ", "_")
        }
        if (!writeFiles.containsKey(cur_id)) {
            writeFiles[cur_id] = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(Files.newOutputStream(Paths.get(prefix + "_" + cur_id + ".fastq.gz")))))
        }
        writeFiles[cur_id].println(line)
    }

    input.close()
    writeFiles.each{k,v -> v.flush(); v.close()}
}
