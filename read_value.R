read_value<-function(prompt_text = "", prompt_suffix = getOption("prompt"), coerce_to="character")
{
  prompt <-paste(prompt_text, prompt_suffix)
  answer<-as(readline(prompt), coerce_to)
  print(answer)
}