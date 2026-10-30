-- Permit long code filenames to wrap within PDF table columns and paragraphs.
function Code(element)
  local escaped = {}
  local special = {
    ['\\'] = '\\textbackslash{}', ['{'] = '\\{', ['}'] = '\\}',
    ['#'] = '\\#', ['$'] = '\\$', ['%'] = '\\%', ['&'] = '\\&',
    ['_'] = '\\_', ['^'] = '\\textasciicircum{}',
    ['~'] = '\\textasciitilde{}'
  }
  for character in element.text:gmatch('.') do
    table.insert(escaped, special[character] or character)
  end
  return pandoc.RawInline('latex', '\\texttt{' ..
    table.concat(escaped, '\\allowbreak{}') .. '}')
end

function Table(element)
  -- Support both the older and current Pandoc table representations.
  local widths = {0.18, 0.28, 0.11, 0.13, 0.12, 0.18}
  local first_header = element.headers and element.headers[1]
  if element.head and #element.head.rows > 0 then
    first_header = element.head.rows[1].cells[1].contents
  end
  if first_header and pandoc.utils.stringify(first_header) == 'Package' then
    widths = {0.21, 0.1233, 0.21, 0.1233, 0.21, 0.1234}
  end
  if element.colspecs and #element.colspecs == 6 then
    for index, width in ipairs(widths) do
      element.colspecs[index][2] = width
    end
  elseif element.widths and #element.widths == 6 then
    element.widths = widths
  end
  return element
end

function Header(element)
  if pandoc.utils.stringify(element.content) == 'Appendix: R Package Versions' then
    return {pandoc.RawBlock('latex', '\\clearpage'), element}
  end
end
