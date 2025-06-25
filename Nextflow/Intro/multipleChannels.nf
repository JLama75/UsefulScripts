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
  x = channel.of(1, 2)
  y = channel.of('a', 'b', 'c')
  echo(x, y)
}

//channel values are consumed sequentially and any empty channel will cause the process to wait, even if the other channels have values.
//The process echo is executed two times because the x channel emits only two values, therefore the c element is discarded. It outputs:
// 1 and a
// 2 and b
