ReadValue <- function(prompt_text = "", prompt_suffix = getOption("prompt"), coerce_to="character")
  # Function reads input from user
  #
  # Args:
  # prompt_text: Text to be shown to user.
  # prompt_suffix: A non-empty string to be used for R's prompt.
  # coerce_to: String specifying the type the input will be coerced to.
  #
{
  prompt <- paste(prompt_text, prompt_suffix)
  answer <- as(readline(prompt), coerce_to)
  print(answer)
}
