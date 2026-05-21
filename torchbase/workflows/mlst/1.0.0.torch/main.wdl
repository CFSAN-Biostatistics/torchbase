version 1.0

workflow mlst {
  input {
    File? reads
    File? contigs
  }

  output {
    String result = "pending"
  }
}
