#' Academic Journal Color Palettes
#'
#' Color palettes inspired by academic journals and scientific publications
#' including Nature, Science, Cell, JACS, and other prestigious journals.
#' Designed specifically for scientific data visualization.
#'
#' @section Available Palettes:
#' \describe{
#'   \item{sci}{General scientific palette with balanced colors for most applications}
#'   \item{nature}{Nature journal inspired colors}
#'   \item{science}{Science journal inspired colors}
#'   \item{cell}{Cell journal inspired colors}
#'   \item{jacs}{JACS journal inspired colors}
#'   \item{fuel}{Fuel journal inspired colors with high contrast}
#'   \item{chem_eng}{Chemical Engineering Journal inspired colors}
#'   \item{nat_comm}{Nature Communications inspired colors}
#'   \item{shinkai}{Makoto Shinkai (新海诚) animation inspired colors}
#' }
#'
#' @name academic_palettes
NULL

# Academic journal color database
academic_db <- list(
  "sci" = list(
    "default" = c("#0C6291", "#A63446", "#F6B101", "#53D600", "#7000AC", "#FF6B6B", 
                  "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F06292", "#5E82A2")
  ),
  "nature" = list(
    "default" = c("#0074B3", "#EC3232", "#F6B101", "#53D600", "#7000AC", "#FFA500", 
                  "#3C5488", "#FB6467", "#526E2D", "#80C4BB", "#ED8474", "#AE5D9D"),
    "materials" = c("#FFA810", "#7000AC", "#53D600", "#D61737", "#F7BF25", "#3634F2",
                    "#5A08A5", "#16386A", "#1E69B0", "#64B0DF", "#770020", "#69B78F")
  ),
  "science" = list(
    "default" = c("#B42B22", "#315A89", "#996E2E", "#EC3232", "#0787C3", "#F6944B",
                  "#8FC0A9", "#467F79", "#C6133B", "#90162D", "#93A5CB", "#2E4F4A")
  ),
  "cell" = list(
    "default" = c("#427AB2", "#F09148", "#FF9896", "#299D8F", "#E9C46A", "#D87659",
                  "#DBDB8D", "#C59D94", "#AFC7E8", "#EA8379", "#7DAEE0", "#B395BD")
  ),
  "jacs" = list(
    "default" = c("#D15354", "#5094D5", "#E8B86C", "#8887CB", "#5E82A2", "#BFC7E5",
                  "#F9AD95", "#ABD8E5", "#F6AD8A", "#3DA6AE", "#F79647", "#C3D69C")
  ),
  "fuel" = list(
    "default" = c("#000035", "#730101", "#009E2B", "#FBC40F", "#C34A00", "#74379F",
                  "#1D4201", "#0172BE", "#DA5319", "#EDB11F")
  ),
  "chem_eng" = list(
    "default" = c("#EB0000", "#0000E4", "#155C00", "#E97200", "#295ABA", "#D47232",
                  "#379124", "#077DF0", "#FF2F4C", "#737373", "#696969", "#242424")
  ),
  "nat_comm" = list(
    "default" = c("#7E4909", "#0E8585", "#830783", "#FA5454", "#C31D1D", "#44988F",
                  "#EBB34F", "#8389C2", "#D87A8A", "#92CEC8", "#ED8666", "#8036FC")
  ),
  "shinkai" = list(
    "default" = c("#5A5FA3", "#8F97C9", "#BFC0DE", "#3F5B72", "#6178A5", "#79A0B5",
                  "#1E3F66", "#2E5C8A", "#4178A5", "#5F8FB8", "#7BA3C7", "#2A4D66"),
    "blue_tones" = c("#1E3F66", "#2E5C8A", "#4178A5", "#5F8FB8", "#7BA3C7", "#2A4D66",
                     "#3C6278", "#4F7A8F", "#6391A5", "#5C7B88", "#7DA2B0", "#9BB8C6"),
    "red_tones" = c("#710000", "#980000", "#C60000", "#E50000", "#F80000")
  )
)

# Continuous palette colors for gradient scales
academic_continuous <- list(
  "sci_gradient" = c("#0C6291", "#FBFEF9", "#A63446"),
  "nature_gradient" = c("#0074B3", "#F0F8FF", "#EC3232"),
  "science_gradient" = c("#315A89", "#F5F5F5", "#B42B22"),
  "cell_gradient" = c("#427AB2", "#FFFFFF", "#F09148"),
  "jacs_gradient" = c("#5094D5", "#FFEEE0", "#D15354"),
  "fuel_gradient" = c("#000035", "#EEEEEE", "#730101"),
  "chem_eng_gradient" = c("#0000E4", "#F0F0F0", "#EB0000"),
  "nat_comm_gradient" = c("#0E8585", "#F8F8F8", "#7E4909"),
  "shinkai_gradient" = c("#1E3F66", "#E0E2F2", "#5A5FA3")
)

# Helper function to check ggplot2 version (from ggsci)
is_ggplot2_350 <- function() {
  utils::packageVersion("ggplot2") >= "3.5.0"
}

#' Academic journal color palettes
#'
#' Color palettes inspired by plots in academic journals.
#' Designed for scientific data visualization with carefully selected colors
#' that work well for both print and digital formats.
#'
#' @param palette Palette type. Available options depend on the journal:
#'   \itemize{
#'     \item sci: "default" (12-color balanced scientific palette)
#'     \item nature: "default", "materials" (12-color Nature-inspired palettes)
#'     \item science: "default" (12-color Science journal palette)
#'     \item cell: "default" (12-color Cell journal palette)
#'     \item jacs: "default" (12-color JACS palette)
#'     \item fuel: "default" (10-color high-contrast palette)
#'     \item chem_eng: "default" (12-color Chemical Engineering Journal palette)
#'     \item nat_comm: "default" (12-color Nature Communications palette)
#'     \item shinkai: "default", "blue_tones", "red_tones" (Makoto Shinkai inspired)
#'   }
#' @param alpha Transparency level, a real number in (0, 1].
#'   See `alpha` in [grDevices::rgb()] for details.
#'
#' @export pal_sci
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom scales manual_pal
#'
#' @examples
#' library("scales")
#' show_col(pal_sci("default")(10))
#' show_col(pal_sci("default", alpha = 0.6)(10))
pal_sci <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"sci"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_nature <- function(palette = c("default", "materials"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"nature"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_science <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"science"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_cell <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"cell"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_jacs <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"jacs"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_fuel <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"fuel"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_chem_eng <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"chem_eng"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_nat_comm <- function(palette = c("default"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"nat_comm"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

#' @rdname pal_sci
#' @export
pal_shinkai <- function(palette = c("default", "blue_tones", "red_tones"), alpha = 1) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  raw_cols <- academic_db$"shinkai"[[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  manual_pal(unname(alpha_cols))
}

# Discrete color scales ======================================================

#' Academic journal color scales
#'
#' See [pal_sci()] and related palette functions for details.
#'
#' @inheritParams pal_sci
#' @param ... Additional parameters for [ggplot2::discrete_scale()].
#'
#' @export scale_color_sci
#'
#' @importFrom ggplot2 discrete_scale
#'
#' @examples
#' library("ggplot2")
#' data("iris")
#'
#' # Basic usage
#' ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#'   geom_point(size = 2) +
#'   scale_color_sci()
#'
#' # With transparency
#' ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#'   geom_point(size = 2) +
#'   scale_color_sci(alpha = 0.7)
scale_color_sci <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_sci(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "sci", palette = pal_sci(palette, alpha), ...)
  }
}

#' @export scale_colour_sci
#' @rdname scale_color_sci
scale_colour_sci <- scale_color_sci

#' @export scale_fill_sci
#' @importFrom ggplot2 discrete_scale
#' @rdname scale_color_sci
scale_fill_sci <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_sci(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "sci", palette = pal_sci(palette, alpha), ...)
  }
}

# Nature scales
#' @export scale_color_nature
#' @rdname scale_color_sci
scale_color_nature <- function(palette = c("default", "materials"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_nature(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "nature", palette = pal_nature(palette, alpha), ...)
  }
}

#' @export scale_colour_nature
#' @rdname scale_color_sci
scale_colour_nature <- scale_color_nature

#' @export scale_fill_nature
#' @rdname scale_color_sci
scale_fill_nature <- function(palette = c("default", "materials"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_nature(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "nature", palette = pal_nature(palette, alpha), ...)
  }
}

# Science scales
#' @export scale_color_science
#' @rdname scale_color_sci
scale_color_science <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_science(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "science", palette = pal_science(palette, alpha), ...)
  }
}

#' @export scale_colour_science
#' @rdname scale_color_sci
scale_colour_science <- scale_color_science

#' @export scale_fill_science
#' @rdname scale_color_sci
scale_fill_science <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_science(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "science", palette = pal_science(palette, alpha), ...)
  }
}

# Cell scales
#' @export scale_color_cell
#' @rdname scale_color_sci
scale_color_cell <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_cell(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "cell", palette = pal_cell(palette, alpha), ...)
  }
}

#' @export scale_colour_cell
#' @rdname scale_color_sci
scale_colour_cell <- scale_color_cell

#' @export scale_fill_cell
#' @rdname scale_color_sci
scale_fill_cell <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_cell(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "cell", palette = pal_cell(palette, alpha), ...)
  }
}

# JACS scales
#' @export scale_color_jacs
#' @rdname scale_color_sci
scale_color_jacs <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_jacs(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "jacs", palette = pal_jacs(palette, alpha), ...)
  }
}

#' @export scale_colour_jacs
#' @rdname scale_color_sci
scale_colour_jacs <- scale_color_jacs

#' @export scale_fill_jacs
#' @rdname scale_color_sci
scale_fill_jacs <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_jacs(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "jacs", palette = pal_jacs(palette, alpha), ...)
  }
}

# Fuel scales
#' @export scale_color_fuel
#' @rdname scale_color_sci
scale_color_fuel <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_fuel(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "fuel", palette = pal_fuel(palette, alpha), ...)
  }
}

#' @export scale_colour_fuel
#' @rdname scale_color_sci
scale_colour_fuel <- scale_color_fuel

#' @export scale_fill_fuel
#' @rdname scale_color_sci
scale_fill_fuel <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_fuel(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "fuel", palette = pal_fuel(palette, alpha), ...)
  }
}

# Chemical Engineering Journal scales
#' @export scale_color_chem_eng
#' @rdname scale_color_sci
scale_color_chem_eng <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_chem_eng(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "chem_eng", palette = pal_chem_eng(palette, alpha), ...)
  }
}

#' @export scale_colour_chem_eng
#' @rdname scale_color_sci
scale_colour_chem_eng <- scale_color_chem_eng

#' @export scale_fill_chem_eng
#' @rdname scale_color_sci
scale_fill_chem_eng <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_chem_eng(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "chem_eng", palette = pal_chem_eng(palette, alpha), ...)
  }
}

# Nature Communications scales
#' @export scale_color_nat_comm
#' @rdname scale_color_sci
scale_color_nat_comm <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_nat_comm(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "nat_comm", palette = pal_nat_comm(palette, alpha), ...)
  }
}

#' @export scale_colour_nat_comm
#' @rdname scale_color_sci
scale_colour_nat_comm <- scale_color_nat_comm

#' @export scale_fill_nat_comm
#' @rdname scale_color_sci
scale_fill_nat_comm <- function(palette = c("default"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_nat_comm(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "nat_comm", palette = pal_nat_comm(palette, alpha), ...)
  }
}

# Shinkai scales
#' @export scale_color_shinkai
#' @rdname scale_color_sci
scale_color_shinkai <- function(palette = c("default", "blue_tones", "red_tones"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("colour", palette = pal_shinkai(palette, alpha), ...)
  } else {
    discrete_scale("colour", scale_name = "shinkai", palette = pal_shinkai(palette, alpha), ...)
  }
}

#' @export scale_colour_shinkai
#' @rdname scale_color_sci
scale_colour_shinkai <- scale_color_shinkai

#' @export scale_fill_shinkai
#' @rdname scale_color_sci
scale_fill_shinkai <- function(palette = c("default", "blue_tones", "red_tones"), alpha = 1, ...) {
  palette <- match.arg(palette)
  if (is_ggplot2_350()) {
    discrete_scale("fill", palette = pal_shinkai(palette, alpha), ...)
  } else {
    discrete_scale("fill", scale_name = "shinkai", palette = pal_shinkai(palette, alpha), ...)
  }
}

# Continuous color scales ====================================================

#' Continuous academic journal color scales
#'
#' Continuous color scales inspired by academic journals, designed for
#' heatmaps, continuous data visualization, and gradient mapping.
#'
#' @param palette Palette type for continuous scales:
#'   \itemize{
#'     \item "sci_gradient": Scientific blue-white-red gradient
#'     \item "nature_gradient": Nature journal blue-white-red gradient
#'     \item "science_gradient": Science journal blue-white-red gradient
#'     \item "cell_gradient": Cell journal blue-white-orange gradient
#'     \item "jacs_gradient": JACS blue-cream-red gradient
#'     \item "fuel_gradient": Fuel journal blue-gray-red gradient
#'     \item "chem_eng_gradient": Chemical Engineering blue-gray-red gradient
#'     \item "nat_comm_gradient": Nature Communications cyan-gray-brown gradient
#'     \item "shinkai_gradient": Shinkai blue-lavender-purple gradient
#'   }
#' @param alpha Transparency level, a real number in (0, 1].
#' @param reverse Logical. Should the order of the colors be reversed?
#' @param ... Additional parameters for [ggplot2::scale_color_gradientn()] or
#'   [ggplot2::scale_fill_gradientn()].
#'
#' @export scale_color_sci_c
#'
#' @importFrom ggplot2 scale_color_gradientn scale_fill_gradientn
#' @importFrom grDevices colorRamp rgb
#'
#' @examples
#' library("ggplot2")
#' 
#' # Create sample data
#' data <- expand.grid(x = 1:10, y = 1:10)
#' data$z <- with(data, x * y)
#'
#' # Continuous color scale
#' ggplot(data, aes(x = x, y = y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_sci_c() +
#'   theme_minimal()
#'
#' # Different palette
#' ggplot(data, aes(x = x, y = y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_sci_c(palette = "nature_gradient") +
#'   theme_minimal()
#'
#' # Science journal style
#' ggplot(data, aes(x = x, y = y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_science_c() +
#'   theme_minimal()
scale_color_sci_c <- function(palette = c("sci_gradient", "nature_gradient", "science_gradient", "cell_gradient", "jacs_gradient", "fuel_gradient", "chem_eng_gradient", "nat_comm_gradient", "shinkai_gradient"), 
                              alpha = 1, reverse = FALSE, ...) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  
  colors <- academic_continuous[[palette]]
  if (reverse) colors <- rev(colors)
  
  # Create smooth gradient
  func_cols <- colorRamp(colors, space = "Lab", interpolate = "spline")
  mat_cols <- func_cols(seq(0L, 1L, length.out = 256))
  gradient_cols <- rgb(
    mat_cols[, 1L], mat_cols[, 2L], mat_cols[, 3L],
    alpha = alpha * 255L, maxColorValue = 255L
  )
  
  scale_color_gradientn(colours = gradient_cols, ...)
}

#' @export scale_colour_sci_c
#' @rdname scale_color_sci_c
scale_colour_sci_c <- scale_color_sci_c

#' @export scale_fill_sci_c
#' @rdname scale_color_sci_c
scale_fill_sci_c <- function(palette = c("sci_gradient", "nature_gradient", "science_gradient", "cell_gradient", "jacs_gradient", "fuel_gradient", "chem_eng_gradient", "nat_comm_gradient", "shinkai_gradient"), 
                             alpha = 1, reverse = FALSE, ...) {
  palette <- match.arg(palette)
  if (alpha > 1L || alpha <= 0L) stop("alpha must be in (0, 1]")
  
  colors <- academic_continuous[[palette]]
  if (reverse) colors <- rev(colors)
  
  # Create smooth gradient
  func_cols <- colorRamp(colors, space = "Lab", interpolate = "spline")
  mat_cols <- func_cols(seq(0L, 1L, length.out = 256))
  gradient_cols <- rgb(
    mat_cols[, 1L], mat_cols[, 2L], mat_cols[, 3L],
    alpha = alpha * 255L, maxColorValue = 255L
  )
  
  scale_fill_gradientn(colours = gradient_cols, ...)
}

# Additional continuous scales for other journals
#' @export scale_color_nature_c
#' @rdname scale_color_sci_c
scale_color_nature_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "nature_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_nature_c
#' @rdname scale_color_sci_c
scale_colour_nature_c <- scale_color_nature_c

#' @export scale_fill_nature_c
#' @rdname scale_color_sci_c
scale_fill_nature_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "nature_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_color_cell_c
#' @rdname scale_color_sci_c
scale_color_cell_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "cell_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_cell_c
#' @rdname scale_color_sci_c
scale_colour_cell_c <- scale_color_cell_c

#' @export scale_fill_cell_c
#' @rdname scale_color_sci_c
scale_fill_cell_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "cell_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_color_jacs_c
#' @rdname scale_color_sci_c
scale_color_jacs_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "jacs_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_jacs_c
#' @rdname scale_color_sci_c
scale_colour_jacs_c <- scale_color_jacs_c

#' @export scale_fill_jacs_c
#' @rdname scale_color_sci_c
scale_fill_jacs_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "jacs_gradient", alpha = alpha, reverse = reverse, ...)
}

# Science continuous scales
#' @export scale_color_science_c
#' @rdname scale_color_sci_c
scale_color_science_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "science_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_science_c
#' @rdname scale_color_sci_c
scale_colour_science_c <- scale_color_science_c

#' @export scale_fill_science_c
#' @rdname scale_color_sci_c
scale_fill_science_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "science_gradient", alpha = alpha, reverse = reverse, ...)
}

# Fuel continuous scales
#' @export scale_color_fuel_c
#' @rdname scale_color_sci_c
scale_color_fuel_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "fuel_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_fuel_c
#' @rdname scale_color_sci_c
scale_colour_fuel_c <- scale_color_fuel_c

#' @export scale_fill_fuel_c
#' @rdname scale_color_sci_c
scale_fill_fuel_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "fuel_gradient", alpha = alpha, reverse = reverse, ...)
}

# Chemical Engineering continuous scales
#' @export scale_color_chem_eng_c
#' @rdname scale_color_sci_c
scale_color_chem_eng_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "chem_eng_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_chem_eng_c
#' @rdname scale_color_sci_c
scale_colour_chem_eng_c <- scale_color_chem_eng_c

#' @export scale_fill_chem_eng_c
#' @rdname scale_color_sci_c
scale_fill_chem_eng_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "chem_eng_gradient", alpha = alpha, reverse = reverse, ...)
}

# Nature Communications continuous scales
#' @export scale_color_nat_comm_c
#' @rdname scale_color_sci_c
scale_color_nat_comm_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "nat_comm_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_nat_comm_c
#' @rdname scale_color_sci_c
scale_colour_nat_comm_c <- scale_color_nat_comm_c

#' @export scale_fill_nat_comm_c
#' @rdname scale_color_sci_c
scale_fill_nat_comm_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "nat_comm_gradient", alpha = alpha, reverse = reverse, ...)
}

# Shinkai continuous scales
#' @export scale_color_shinkai_c
#' @rdname scale_color_sci_c
scale_color_shinkai_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_color_sci_c(palette = "shinkai_gradient", alpha = alpha, reverse = reverse, ...)
}

#' @export scale_colour_shinkai_c
#' @rdname scale_color_sci_c
scale_colour_shinkai_c <- scale_color_shinkai_c

#' @export scale_fill_shinkai_c
#' @rdname scale_color_sci_c
scale_fill_shinkai_c <- function(alpha = 1, reverse = FALSE, ...) {
  scale_fill_sci_c(palette = "shinkai_gradient", alpha = alpha, reverse = reverse, ...)
}