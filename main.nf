params.str = 'Hello world!'

process splitLetters {
  container 'nextflow/nextflow:21.10.6'
  output:
    path 'chunk_*'

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

process convertToUpper {
  container 'nextflow/nextflow:21.10.6'
  input:
    path x
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}

workflow {
  splitLetters | flatten | convertToUpper | view { it.trim() }
}

