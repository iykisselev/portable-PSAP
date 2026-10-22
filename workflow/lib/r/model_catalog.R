# Disease models and status payloads come from the same catalog as generation.

catalog_strings <- function(value, label, nonempty = TRUE) {
  valid <- is.list(value) && all(vapply(value, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1)))
  if (!valid || (nonempty && !length(value)) || anyDuplicated(unlist(value))) {
    stop("Catalog ", label, " must be an array of unique nonempty strings", call. = FALSE)
  }
  as.character(unlist(value, use.names = FALSE))
}


read_model_catalog <- function(path) {
  catalog <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!is.list(catalog) || !identical(catalog$schema_version, 1L)) stop("Unsupported model catalog schema")
  catalog$roles <- catalog_strings(catalog$roles, "roles")
  if (!is.list(catalog$models) || !length(catalog$models)) stop("Catalog models must be a nonempty array")
  identifiers <- tags <- character()
  for (i in seq_along(catalog$models)) {
    model <- catalog$models[[i]]
    if (!is.list(model)) stop("Catalog model entries must be objects")
    for (name in c("id", "role", "scope", "genotype_class", "tag")) {
      value <- model[[name]]
      if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
        stop("Catalog model requires ", name)
      }
    }
    if (model$id %in% identifiers || model$tag %in% tags || !model$role %in% catalog$roles ||
        !model$scope %in% c("autosomal", "non_PAR_chrX") || !is.logical(model$requires_hemizygote) ||
        length(model$requires_hemizygote) != 1L || is.na(model$requires_hemizygote)) {
      stop("Invalid or duplicate catalog model ", model$id)
    }
    identifiers <- c(identifiers, model$id)
    tags <- c(tags, model$tag)
    model$compatible_control_models <- catalog_strings(model$compatible_control_models, "compatible_control_models")
    catalog$models[[i]] <- model
  }
  if (any(vapply(catalog$models, function(m) any(!m$compatible_control_models %in% identifiers), logical(1)))) {
    stop("Catalog contains an unknown compatible control model")
  }
  if (!is.list(catalog$statuses) || !length(catalog$statuses)) stop("Catalog statuses must be a nonempty array")
  statuses <- list()
  for (status in catalog$statuses) {
    if (!is.list(status)) stop("Catalog status entries must be objects")
    code <- status$code
    if (!is.character(code) || length(code) != 1L || is.na(code) || !nzchar(code) || code %in% names(statuses)) {
      stop("Invalid or duplicate catalog status")
    }
    status$domains <- catalog_strings(status$domains, "status domains")
    if (any(!status$domains %in% c("scope", "role", "partition", "aggregation", "triage"))) {
      stop("Unknown catalog status domain for ", code)
    }
    if (!is.list(status$required_payload) || !setequal(names(status$required_payload), status$domains)) {
      stop("Catalog status ", code, " must declare payloads for every domain")
    }
    status$required_payload <- lapply(
      status$required_payload, catalog_strings,
      label = "required_payload", nonempty = FALSE
    )
    if (any(!unlist(status$required_payload, use.names = FALSE) %in% c("role", "failed_roles"))) {
      stop("Unknown catalog payload field for ", code)
    }
    if (!is.logical(status$diagnostic_only) || length(status$diagnostic_only) != 1L || is.na(status$diagnostic_only)) {
      stop("Catalog status ", code, " requires diagnostic_only boolean")
    }
    statuses[[code]] <- status
  }
  catalog$statuses <- statuses
  catalog
}


enabled_catalog_models <- function(hemizygote) {
  if (!is.logical(hemizygote) || length(hemizygote) != 1L || is.na(hemizygote)) stop("hemizygote must be a boolean")
  Filter(function(model) hemizygote || !model$requires_hemizygote, PSAP_MODEL_CATALOG$models)
}


model_roles <- function(hemizygote) {
  used <- vapply(enabled_catalog_models(hemizygote), function(model) model$role, character(1))
  PSAP_MODEL_CATALOG$roles[PSAP_MODEL_CATALOG$roles %in% used]
}


validate_status <- function(code, domain, role = NULL, failed_roles = character(), production = TRUE) {
  if (!is.character(code) || length(code) != 1L || is.na(code) ||
      !is.character(domain) || length(domain) != 1L || is.na(domain)) stop("Invalid status code or domain")
  definition <- PSAP_MODEL_CATALOG$statuses[[code]]
  if (is.null(definition) || !domain %in% definition$domains) stop("Unsupported ", domain, " status '", code, "'")
  if (production && definition$diagnostic_only) stop("Status '", code, "' is diagnostic only")
  if (!is.null(role) && (!is.character(role) || length(role) != 1L || is.na(role) ||
      !role %in% PSAP_MODEL_CATALOG$roles)) stop("Unknown diagnostic role")
  if (!is.character(failed_roles) || anyNA(failed_roles) || anyDuplicated(failed_roles) ||
      !identical(failed_roles, PSAP_MODEL_CATALOG$roles[PSAP_MODEL_CATALOG$roles %in% failed_roles])) {
    stop("failed_roles must contain unique roles in catalog order")
  }
  required <- definition$required_payload[[domain]]
  if ("role" %in% required && is.null(role)) stop("Status '", code, "' requires role")
  if ("failed_roles" %in% required && !length(failed_roles)) stop("Status '", code, "' requires nonempty failed_roles")
  if (length(failed_roles) && !"failed_roles" %in% required) {
    stop("Status '", code, "' does not allow failed_roles in ", domain)
  }
  if (!is.null(role) && !"role" %in% required) stop("Status '", code, "' does not allow role in ", domain)
  invisible(code)
}
