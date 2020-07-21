#!/usr/bin/env groovy

@Grapes(
  @Grab(group='io.reactivex.rxjava3', module='rxjava', version='3.0.4')
)
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.nio.file.Files
import java.nio.file.Paths
import java.util.regex.Matcher
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
cli.r(longOpt:'renumber', args:0, argName: 'renumber', 'Renumber the reads')
cli.h(longOpt:'help', 'Show help')

def options = cli.parse(args)

if (!options || options.h || !options.f || !options.p) {
  println(cli.usage())
  System.exit(1)
}

String file = options.f
String prefix = options.p
process(file, prefix, options.r)

class Entry {
    String id;
    BufferedWriter writer;
    int cur_line;
}

@CompileStatic
void process(String file, String prefix, boolean renumber) {
    BufferedReader input = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(Paths.get(file)))))

    Map<String, Entry> writeFiles = [:]

    String cur_id = null
    input.eachLine{line, number ->
        if ((number-1) % 4 == 0) {
            Matcher m = line =~ /(.+)[ _](\d+)\/(\d+)$/
            m.find()
            String id = m.group(1)
            String rn = m.group(2)
            String pn = m.group(3)
            cur_id = pn
            if (!writeFiles.containsKey(cur_id)) {
                def writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(Files.newOutputStream(Paths.get(prefix + "_" + cur_id + ".fastq.gz")))))
                writeFiles[cur_id] = new Entry('id': cur_id, 'writer': writer, cur_line:1)
            } else {
                writeFiles[cur_id].cur_line = writeFiles[cur_id].cur_line + 1
            }
            if (renumber) {
                line = "${id}_${writeFiles[cur_id].cur_line}/${pn}"
            } else {
                line = "${id}_${rn}/${pn}"
            }
        }
        writeFiles[cur_id].writer.println(line)
    }

    input.close()
    writeFiles.each{k,v -> v.writer.flush(); v.writer.close()}
}
