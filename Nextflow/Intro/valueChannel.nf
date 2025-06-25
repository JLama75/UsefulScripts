process echo {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}

workflow {
  x = channel.value(1)
  y = channel.of('a', 'b', 'c')
  echo(x, y)
}

// By definition, a value channel is bound to a single value and it can be read an unlimited number of times without consuming its content. 
//Therefore, when mixing a value channel with one or more (queue) channels, it does not affect the process termination because the underlying value is applied repeatedly.
//The process termination is determined by the contents of y. It outputs:
//1 and a
//1 and b
//1 and c
